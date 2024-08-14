import os
import subprocess
import unittest

import netCDF4
import numpy as np


class SingleColumnThermo(unittest.TestCase):
    # A few useful global variables for the class
    executable = "../nextsim"

    init_file = "init_column.nc"
    config_file = "ThermoIntegration.cfg"
    diagnostics_file = "diagnostic.nc"

    # Load the data once and re-use for all tests
    hice = np.array([])
    hsnow = np.array([])
    tice = np.array([])

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
        cls.hice = np.squeeze(np.array(root.groups["data"].variables["hice"][:].data))
        cls.hsnow = np.squeeze(np.array(root.groups["data"].variables["hsnow"][:].data))
        cls.tice = np.array(root.groups["data"].variables["tice"][:].data)

    @classmethod
    def __make_cfg_file(cls):
        cfg = open(cls.config_file, "w")
        cfg.write("""
[model]
init_file = init_column.nc
start = 1900-01-01T00:00:00Z
stop = 2011-01-01T00:00:00Z
time_step = P0-1T00:00:00

[Modules]
DiagnosticOutputModule = Nextsim::ConfigOutput
IceAlbedoModule = Nextsim::MU71Albedo
AtmosphereBoundaryModule = Nextsim::MU71Atmosphere
OceanBoundaryModule = Nextsim::FluxConfiguredOcean
IceThermodynamicsModule = Nextsim::ThermoWinton

[ConfigOutput]
start = 2010-01-01T00:00:00Z
field_names = hsnow,hice,tice

[FluxConfiguredOcean]
qio = 2
sss = 35
sst = -1.89

[nextsim_thermo]
I_0 = 0.3
ks = 0.31
        """)
        cfg.close()

    @classmethod
    def __make_init_column(cls):
        import sys
        sys.path.append('../python')

        from make_init_base import initMaker

        nfirst = 1
        nsecond = 1
        nLayers = 3
        resolution = 50

        initializer = initMaker("init_column.nc", nfirst, nsecond, nLayers, resolution, checkZeros=False)

        ice_salinity = 5  # should match Ice::s in constants.hpp
        mu: float = -0.055  # should match Water::mu in constants.hpp
        ocean_temperature = -1.54
        ocean_salinity = ocean_temperature / mu

        initializer.mask[:, :] = 1
        initializer.cice[:, :] = 1.
        initializer.hice[:, :] = 3.
        initializer.hsnow[:, :] = 0.3
        initializer.sss[:, :] = ocean_salinity
        initializer.sst[:, :] = ocean_temperature
        initializer.tice[:, :, :] = ice_salinity * mu

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

        mean = 3.1189
        max = 3.3419
        min = 2.9805
        self.assertAlmostEqual(max, self.hice.max(), 4, "Max ice thickness not ~= " + str(max) + " m")
        self.assertAlmostEqual(min, self.hice.min(), 4, "Min ice thickness not ~= " + str(min) + " m")
        self.assertAlmostEqual(mean, self.hice.mean(), 4, "Mean ice thickness not ~= " + str(mean) + " m")

    def test_snowThickness(self):
        """
        Test the snow thickness against standard max, min, and mean values
        """

        mean = 0.2474
        max = 0.4000
        min = 0.0000
        self.assertAlmostEqual(max, self.hsnow.max(), 4, "Max snow thickness not ~= " + str(max) + " m")
        self.assertAlmostEqual(min, self.hsnow.min(), 4, "Min snow thickness not ~= " + str(min) + " m")
        self.assertAlmostEqual(mean, self.hsnow.mean(), 4, "Mean snow thickness not ~= " + str(mean) + " m")

    def test_temperatureTest(self):
        """
        Test the surface and internal temperatures against standard max, min, and mean values

        NB! Here, I put the "places" argument of assertAlmostEqual to 3 for the mean and min comparison. I do this
        because I get inconsistent results on different platforms in the GitHub CI(!) The reason is that the testing
        framework compares up to a given number of decimal places, but the mean and min temperatures have two
        significant digits before the decimal point, making the comparison for those equal to comparing six
        significant digits, when the assertAlmostEqual "places" argument is set to 4. In the GitHub CI the seventh
        significant digit changes between 4 and 5 for the T1 mean, so the result is either -17.6250 or -17.6249 - up
        to 4 digits. This is normal, because the output is only accurate to six significant digits anyway.
        """

        mean = [-17.6250, -7.6068, -3.7998]
        max = [0.0000, -1.1336, -1.5975]
        min = [-33.1612, -14.8637, -6.1424]
        for i in range(3):
            self.assertAlmostEqual(max[i], self.tice[:, i].max(), 4, "Max T" + str(i) + " not ~= " + str(max[i]) + " m")
            self.assertAlmostEqual(min[i], self.tice[:, i].min(), 3, "Min T" + str(i) + " not ~= " + str(min[i]) + " m")
            self.assertAlmostEqual(mean[i], self.tice[:, i].mean(), 3,
                                   "Mean T" + str(i) + " not ~= " + str(mean[i]) + " m")


if __name__ == '__main__':
    unittest.main()
