import os
import subprocess
import sys
import unittest

import netCDF4
import numpy as np


class Polynya(unittest.TestCase):
    # A few useful global variables for the class
    executable = "../nextsim"

    init_file = "init_polynya.nc"
    config_file = "PolynyaIntegration.cfg"
    diagnostics_file = "polynya.diagnostic.nc"

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
        cls.__make_init_column()

        # Create the config file
        cls.__make_cfg_file()

        # Run the model
        subprocess.run(cls.executable + " --config-file " + cls.config_file, shell=True, check=True)

        # Load the basic variables
        root = netCDF4.Dataset(cls.diagnostics_file, "r", format="NETCDF4")
        cls.cice = np.squeeze(np.array(root.groups["data"].variables["cice"][:].data))
        cls.hice = np.squeeze(np.array(root.groups["data"].variables["hice"][:].data))
        cls.uice = np.array(root.groups["data"].variables["u"][:].data)
        cls.vice = np.array(root.groups["data"].variables["v"][:].data)

    @classmethod
    def __make_cfg_file(cls):
        cfg = open(cls.config_file, "w")
        cfg.write("""
[model]
init_file = init_polynya.nc
start = 2023-01-01T00:00:00Z
stop = 2023-01-05T00:00:00Z
time_step = P0-0T00:30:00
missing_value = 1e20

[Modules]
DiagnosticOutputModule = Nextsim::ConfigOutput
DynamicsModule = Nextsim::BBMDynamics
AtmosphereBoundaryModule = Nextsim::FluxConfiguredAtmosphere
OceanBoundaryModule = Nextsim::FluxConfiguredOcean
IceThermodynamicsModule = Nextsim::ThermoWinton

[ConfigOutput]
period = P0-0T00:01:00
field_names = hice,cice,u,v
filename = polynya.diagnostic.nc

[FluxConfiguredOcean]
qio = 2.30
sss = 28
sst = -1.54
mld = 10.
current_u = 0.
current_v = 0.

[FluxConfiguredAtmosphere]
Q_ia = 30
Q_ow = 300
dQia_dT = 0
wind_u = 16
wind_v = 12
        """)
        cfg.close()

    @classmethod
    def __make_init_column(cls):
        sys.path.append('../run')
        sys.path.append('../../run')
        from make_init_base import initMaker

        # Creates initial conditions for the Bjornsson et al. (2001) polynya case

        # Domain size [km]
        x = 100
        y = 50
        res = 4


        nfirst = int(y / res)
        nsecond = int(x / res)

        fname = f"init_polynya.nc"

        # The model expects everything in metres
        initializer = initMaker(fname, nfirst, nsecond, res*1e3, isWinton=True, checkZeros=False)

        # Ice everywhere and all boundaries closed, except the x = 100 km end
        initializer.mask[:, :] = 1.
        initializer.mask[0, :] = 0.
        initializer.mask[-1, :] = 0.
        initializer.mask[:, 0] = 0.
        #initializer.mask[:, -1] = 0. ## right

        # Uniform concentration of 90%
        initializer.cice[:, :] = 0.9

        # Uniform thickness of 20 cm
        initializer.hice[:, :] = 0.2

        # Undamaged ice
        initializer.damage[:, :] = 1.

        # Ice and ocean temperature and salinity at the freezing point
        ice_salinity = 5  # should match Ice::s in constants.hpp
        mu: float = -0.055  # should match Water::mu in constants.hpp
        ocean_temperature = -1.54
        ocean_salinity = ocean_temperature / mu

        initializer.sss[:, :] = ocean_salinity
        initializer.sst[:, :] = ocean_temperature
        initializer.tsurf[:, :] = ice_salinity * mu
        initializer.tbott = initializer.tsurf
        initializer.tintr = initializer.tsurf

        # All other variables are zero or not needed

        # The file is written when initializer goes out of scope

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

        mean = 0.0301
        max = 2.3082
        min = -1.7963
        self.assertAlmostEqual(max, self.hice.max(), 4, "Max ice thickness not ~= " + str(max) + " m (full DG field)")
        self.assertAlmostEqual(min, self.hice.min(), 4, "Min ice thickness not ~= " + str(min) + " m (full DG field)")
        self.assertAlmostEqual(mean, self.hice.mean(), 4, "Mean ice thickness not ~= " + str(mean) + " m (full DG field)")

    def test_concentration(self):
        """
        Test the ice concentration against standard max, min, and mean values
        """

        mean = 0.0931
        max = 2.2597
        min = -1.6293
        self.assertAlmostEqual(max, self.cice.max(), 4, "Max conentration not ~= " + str(max) + " (full DG field)")
        self.assertAlmostEqual(min, self.cice.min(), 4, "Min concentration not ~= " + str(min) + " (full DG field)")
        self.assertAlmostEqual(mean, self.cice.mean(), 4, "Mean concentration not ~= " + str(mean) + " (full DG field)")

    def test_veloities(self):
        """
        Test the ice velocities against standard max, min, and mean values
        """

        mean = [0.2237, 0.0548]
        max = [0.3910, 0.1305]
        min = [-0.0700, -0.1559]
        for (var, i, name) in zip([self.uice, self.vice], [0,1], ["uice", "vice"]):
            self.assertAlmostEqual(max[i], var.max(), 4, "Max "+name+" velocity not ~= " + str(max[i]) + " m/s (full DG field)")
            self.assertAlmostEqual(min[i], var.min(), 4, "Min "+name+" velocity not ~= " + str(min[i]) + " m/s (full DG field)")
            self.assertAlmostEqual(mean[i], var.mean(), 4, "Mean "+name+" velocity not ~= " + str(mean[i]) + " m/s (full DG field)")


if __name__ == '__main__':
    unittest.main()
