import os
import subprocess
import time
import unittest

import netCDF4
import numpy as np


class SingleColumnThermo(unittest.TestCase):
    """A test case class for a single column thermodynamics model."""

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
        """Set up the test case before running a test."""
        # Make the init column
        cls.__make_init_column()

        # Create the config file
        cls.__make_cfg_file()

        # Run the model
        subprocess.run(
            cls.executable + " --config-file " + cls.config_file, shell=True, check=True
        )

        # Load the basic variables
        ncFile = netCDF4.Dataset(cls.diagnostics_file, "r", format="NETCDF4")
        cls.hice = np.squeeze(np.array(ncFile.variables["hice"][:].data))
        cls.hsnow = np.squeeze(np.array(ncFile.variables["hsnow"][:].data))
        cls.tsurf = np.array(ncFile.variables["tsurf"][:].data)
        cls.tintr = np.array(ncFile.variables["tinterior"][:].data)
        cls.tbott = np.array(ncFile.variables["tbottom"][:].data)

    @classmethod
    def __make_cfg_file(cls):
        with open(cls.config_file, "w") as cfg:
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
field_names = hsnow,hice,tsurf,tinterior,tbottom

[FluxConfiguredOcean]
qio = 2
sss = 35
sst = -1.89

[nextsim_thermo]
I_0 = 0.3
ks = 0.31
        """)

    @classmethod
    def __make_init_column(cls):
        ncFile = netCDF4.Dataset(cls.init_file, "w", format="NETCDF4")
        structure_name = "parametric_rectangular"
        ncFile.structure_name = structure_name

        time_var = ncFile.createVariable("time_meta", "i8")
        data_time = 1263204000
        time_var[:] = data_time
        time.units = "seconds since 1970-01-01T00:00:00Z"
        formatted = ncFile.createVariable("formatted", str)
        formatted.format = "%Y-%m-%dT%H:%M:%SZ"
        formatted[0] = "2010-01-01T00:00:00Z"

        nfirst = 1
        nsecond = 1
        ncg = 1
        n_dg = 1
        n_dgstress = 1
        n_coords = 2

        yDim = ncFile.createDimension("ydim", nfirst)
        xDim = ncFile.createDimension("xdim", nsecond)
        yVertexDim = ncFile.createDimension("yvertex", nfirst + 1)
        xVertexDim = ncFile.createDimension("xvertex", nsecond + 1)
        ycg_dim = ncFile.createDimension("y_cg", nfirst * ncg + 1)
        xcg_dim = ncFile.createDimension("x_cg", nsecond * ncg + 1)
        dg_comp = ncFile.createDimension("dg_comp", n_dg)
        dgs_comp = ncFile.createDimension("dgstress_comp", n_dgstress)
        n_coords_comp = ncFile.createDimension("ncoords", n_coords)

        field_dims = ("ydim", "xdim")
        coord_dims = ("yvertex", "xvertex", "ncoords")

        # Array coordinates
        node_lon = np.zeros((nfirst + 1, nsecond + 1))
        node_lat = np.zeros((nfirst + 1, nsecond + 1))

        node_lon[0, 0] = 0
        node_lon[0, 1] = 270
        node_lon[1, 1] = 180
        node_lon[1, 0] = 90

        lat0 = 89.84  # Gives a box around the pole 25 km a side

        node_lat[:, :] = lat0

        coords = ncFile.createVariable("coords", "f8", coord_dims)
        coords[:, :, 0] = node_lon
        coords[:, :, 1] = node_lat

        elem_lon = ncFile.createVariable("longitude", "f8", field_dims)
        elem_lon[:, :] = 0
        elem_lat = ncFile.createVariable("latitude", "f8", field_dims)
        elem_lat[:, :] = 90

        grid_azimuth = ncFile.createVariable("grid_azimuth", "f8", field_dims)
        grid_azimuth[:, :] = 0

        mu: float = -0.055  # should match Water::mu in constants.hpp
        ocean_temperature = -1.54
        ocean_salinity = ocean_temperature / mu

        mask = ncFile.createVariable("mask", "f8", field_dims)
        mask[:, :] = [[1]]
        cice = ncFile.createVariable("cice", "f8", field_dims)
        cice[:, :] = 1.0
        hice = ncFile.createVariable("hice", "f8", field_dims)
        hice[:, :] = 3.00
        hsnow = ncFile.createVariable("hsnow", "f8", field_dims)
        hsnow[:, :] = 0.3
        sss = ncFile.createVariable("sss", "f8", field_dims)
        sss[:, :] = ocean_salinity
        sst = ncFile.createVariable("sst", "f8", field_dims)
        sst[:, :] = ocean_temperature
        # Ice is at rest
        u = ncFile.createVariable("u", "f8", field_dims)
        u[:, :] = 0
        v = ncFile.createVariable("v", "f8", field_dims)
        v[:, :] = 0
        ncFile.close()

    @classmethod
    def tearDownClass(cls):
        """Delete the netCDF output and temporary files."""
        if os.path.isfile(cls.diagnostics_file):
            os.remove(cls.diagnostics_file)

        if os.path.isfile(cls.init_file):
            os.remove(cls.init_file)

        if os.path.isfile(cls.config_file):
            os.remove(cls.config_file)

    def test_iceThickness(self):
        """Test the ice thickness against standard max, min, and mean values."""
        meanval = 3.1093
        maxval = 3.3327
        minval = 2.9702
        hiceDG0 = self.hice[:, 0]
        self.assertAlmostEqual(
            maxval, hiceDG0.max(), 4, f"Max ice thickness not ~= {maxval} m"
        )
        self.assertAlmostEqual(
            minval, hiceDG0.min(), 4, f"Min ice thickness not ~= {minval} m"
        )
        self.assertAlmostEqual(
            meanval, hiceDG0.mean(), 4, f"Mean ice thickness not ~= {meanval} m"
        )

    def test_snowThickness(self):
        """Test the snow thickness against standard max, min, and mean values."""
        meanval = 0.2474
        maxval = 0.4000
        minval = 0.0000
        snowDG0 = self.hsnow[:, 0]
        self.assertAlmostEqual(
            maxval, snowDG0.max(), 4, f"Max snow thickness not ~= {maxval} m"
        )
        self.assertAlmostEqual(
            minval, snowDG0.min(), 4, f"Min snow thickness not ~= {minval} m"
        )
        self.assertAlmostEqual(
            meanval, snowDG0.mean(), 4, f"Mean snow thickness not ~= {meanval} m"
        )

    def test_temperatureTest(self):
        """
        Test the surface and internal temperatures against standard max, min, and mean values.

        NB! Here, I put the "places" argument of assertAlmostEqual to 3 for the mean and min comparison. I do this
        because I get inconsistent results on different platforms in the GitHub CI(!) The reason is that the testing
        framework compares up to a given number of decimal places, but the mean and min temperatures have two
        significant digits before the decimal point, making the comparison for those equal to comparing six
        significant digits, when the assertAlmostEqual "places" argument is set to 4. In the GitHub CI the seventh
        significant digit changes between 4 and 5 for the T1 mean, so the result is either -17.6250 or -17.6249 - up
        to 4 digits. This is normal, because the output is only accurate to six significant digits anyway.
        """
        meanvals = [-17.6202, -7.5904, -3.7944]
        maxvals = [0.0000, -1.1280, -1.5939]
        minvals = [-33.1569, -14.8520, -6.1389]
        for i, t_level in enumerate(
            (self.tsurf[:, 0, 0, 0], self.tintr[:, 0, 0, 0], self.tbott[:, 0, 0, 0])
        ):
            self.assertAlmostEqual(
                maxvals[i], t_level.max(), 4, f"Max T {i} not ~= {maxvals[i]} ˚C"
            )
            self.assertAlmostEqual(
                minvals[i], t_level.min(), 3, f"Min T {i} not ~= {minvals[i]} ˚C"
            )
            self.assertAlmostEqual(
                meanvals[i], t_level.mean(), 3, f"Mean T {i} not ~= {meanvals[i]} ˚C"
            )


if __name__ == "__main__":
    unittest.main()
