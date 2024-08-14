import os
import subprocess
import unittest

import netCDF4
import numpy as np


class Benchmark(unittest.TestCase):
    # A few useful global variables for the class
    executable = "../nextsim"

    init_file = "init_benchmark.nc"
    config_file = "BenchmarkIntegration.cfg"
    diagnostics_file = "benchmark.diagnostic.nc"

    # Load the data once and re-use for all tests
    hice = np.array([])
    cice = np.array([])
    uice = np.array([])
    vice = np.array([])

    @classmethod
    def setUpClass(cls):
        """
        A set-up class which,
          - Creates the initialisation file, using make_init_column.py
          - Runs the model
          - Loads the neccesary variables from the output file
        """

        # Make the init column
        cls.__make_init()

        # Create the config file
        cls.__make_cfg_file()

        # Run the model
        subprocess.run(cls.executable + " --config-file " + cls.config_file, shell=True, check=True)

        # Load the basic variables
        root = netCDF4.Dataset(cls.diagnostics_file, "r", format="NETCDF4")
        cls.hice = np.squeeze(np.array(root.groups["data"].variables["hice"][:].data))[:, 2:-1, 2:-1]
        cls.cice = np.squeeze(np.array(root.groups["data"].variables["cice"][:].data))[:, 2:-1, 2:-1]
        cls.uice = np.squeeze(np.array(root.groups["data"].variables["u"][:].data))[:, 2:-1, 2:-1]
        cls.vice = np.squeeze(np.array(root.groups["data"].variables["v"][:].data))[:, 2:-1, 2:-1]

    @classmethod
    def __make_cfg_file(cls):
        cfg = open(cls.config_file, "w")
        cfg.write("""
[model]
init_file = init_benchmark.nc
start = 2023-01-01T00:00:00Z
stop = 2023-01-03T00:00:00Z
time_step = P0-0T00:02:00

[Modules]
DiagnosticOutputModule = Nextsim::ConfigOutput
DynamicsModule = Nextsim::MEVPDynamics
IceThermodynamicsModule = Nextsim::DummyIceThermodynamics
LateralIceSpreadModule = Nextsim::DummyIceSpread
AtmosphereBoundaryModule = Nextsim::BenchmarkAtmosphere
OceanBoundaryModule = Nextsim::BenchmarkOcean

[ConfigOutput]
period = P0-0T04:00:00
start = 2010-01-01T12:00:00Z
field_names = hsnow,hice,cice,u,v
filename = benchmark.diagnostic.nc
        """)
        cfg.close()

    @classmethod
    def __make_init(cls):
        import sys
        sys.path.append('../python')

        from make_init_base import initMaker
        from math import sin

        L = 512
        res = 16

        nfirst = int(L / res)
        nsecond = int(L / res)
        nLayers = 1

        initializer = initMaker(cls.init_file, nfirst, nsecond, nLayers, res*1e3, checkZeros=False)
        # The model expects everything in metres, while the benchmark problem in Mehlman et al. (2021) is defined in km.

        # Ice everywhere and all boundaries closed
        initializer.mask[:, :] = 1.
        initializer.mask[0, :] = 0.
        initializer.mask[-1, :] = 0.
        initializer.mask[:, 0] = 0.
        initializer.mask[:, -1] = 0.

        # Uniform concentration of 100%
        initializer.cice[:, :] = 1.

        # Loop over ice thickness to construct the initial conditions. This should be a pattern of undulating ice.
        for ix in range(nfirst):
            x = ix * res
            for iy in range(nsecond):
                y = iy * res
                initializer.hice[ix, iy] = 0.3 + 0.005 * (sin(60e-3 * x) + sin(30e-3 * y))

    @classmethod
    def tearDownClass(cls):
        """
        A tear-down class that deletes the netCDF output and temporary files
        """

        if os.path.isfile(cls.diagnostics_file):
            os.remove(cls.diagnostics_file)

        if os.path.isfile(cls.init_file):
            os.remove(cls.init_file)

        if os.path.isfile(cls.config_file):
            os.remove(cls.config_file)

    def test_iceThickness(self):
        """
        Test the ice thickness against standard max, min, and mean values
        """

        mean = 0.2989
        max = 0.3654
        min = 0.1118
        self.assertAlmostEqual(max, self.hice.max(), 4, "Max ice thickness not ~= " + str(max) + " m")
        self.assertAlmostEqual(min, self.hice.min(), 4, "Min ice thickness not ~= " + str(min) + " m")
        self.assertAlmostEqual(mean, self.hice.mean(), 4, "Mean ice thickness not ~= " + str(mean) + " m")

    def test_concentration(self):
        """
        Test the ice concentration against standard max, min, and mean values
        """

        mean = 0.9930
        max = 1.0
        min = 0.6176
        self.assertAlmostEqual(max, self.cice.max(), 4, "Max concentration not ~= " + str(max) + " m")
        self.assertAlmostEqual(min, self.cice.min(), 4, "Min concentration not ~= " + str(min) + " m")
        self.assertAlmostEqual(mean, self.cice.mean(), 4, "Mean concentration not ~= " + str(mean) + " m")

    def test_velocity(self):
        """
        Test the ice velocity against standard max, min, and mean values
        We need a lot of "places" for the minimum velocity
        """

        vel = np.hypot(self.uice,self.vice)

        mean =0.09968
        max = 0.1872
        min = 5.4789e-6
        self.assertAlmostEqual(max, vel.max(), 4, "Max velocity not ~= " + str(max) + " m")
        self.assertAlmostEqual(min, vel.min(), 10, "Min velocity not ~= " + str(min) + " m")
        self.assertAlmostEqual(mean, vel.mean(), 5, "Mean velocity not ~= " + str(mean) + " m")

if __name__ == '__main__':
    unittest.main()
