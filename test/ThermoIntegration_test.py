import os
import subprocess
import unittest

import netCDF4
import numpy as np


class SingleColumnThermo(unittest.TestCase):
    # A few useful global variables for the class
    make_init = "make_init_column.py"
    executable = "./nextsim"

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
        init_script = open(cls.make_init)
        exec(init_script.read())
        init_script.close()

        # Run the model
        subprocess.run(cls.executable + " --config-file " + cls.config_file, shell=True, check=True)

        # Load the basic variables
        root = netCDF4.Dataset(cls.diagnostics_file, "r", format="NETCDF4")
        cls.hice = np.squeeze(np.array(root.groups["data"].variables["hice"][:].data))
        cls.hsnow = np.squeeze(np.array(root.groups["data"].variables["hsnow"][:].data))
        cls.tice = np.array(root.groups["data"].variables["tice"][:].data)

    @classmethod
    def tearDownClass(cls):
        """
        A tear-down class that deletes the netCDF output
        """

        if os.path.isfile(cls.diagnostics_file):
            os.remove(cls.diagnostics_file)

        if os.path.isfile(cls.init_file):
            os.remove(cls.init_file)

    def test_iceThickness(self):
        """
        Test the ice thickness against standard max, min, and mean values
        """

        mean = 3.1326
        max = 3.3554
        min = 2.9947
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
        """

        mean = [-17.6252, -7.5953, -3.7728]
        max = [0.0000, -1.1218, -1.5717]
        min = [-33.1642, -14.8618, -6.1176]
        for i in range(3):
            self.assertAlmostEqual(max[i], self.tice[:, i].max(), 4, "Max T" + str(i) + " not ~= " + str(max[i]) + " m")
            self.assertAlmostEqual(min[i], self.tice[:, i].min(), 4, "Min T" + str(i) + " not ~= " + str(min[i]) + " m")
            self.assertAlmostEqual(mean[i], self.tice[:, i].mean(), 4,
                                   "Mean T" + str(i) + " not ~= " + str(mean[i]) + " m")


if __name__ == '__main__':
    unittest.main()
