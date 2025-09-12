import netCDF4
import numpy as np


class initMaker:
    """
    A "plug-and-play" initialisation class for neXtSIM. The user needs to supply
    at minimum the grid dimensions and resolution. They may also supply any
    initialisation fields they need, as well as a land mask.
    Usage:
     0. Import make_init_base
      >>> from make_init_base import initMaker
     1. Create an initialiser object with a given filename, e.g.
      >>> init = initMaker("test")
         2a. Create a cartesian grid from dimensions and resolution, e.g.
          >>> init.make_cartesian_grid(128, 128, 3e3)
         2b. Create a geographic grid from a file containing coordinates. We need to define type of grid:
            * "p": only center points of the grid are given
            * "ur": the centre and upper right corner of the grid is given
            * "ll": the centre and lower left corner of the grid is given
          >>> init.make_geographic_grid("coords.nc", "p", p_lat_name="latitude", p_lon_name="longitude")
     3. Modify any variables needed, e.g.
      >>> init.cice = 1
      >>> init.hice = 3
      >>> init.hice[10:20, 10:20] = 4
     4. The init file is written when the initialiser object goes out of scope
        (e.g. when the program ends or at the end of a loop).
    """

    def __init__(self, fname, nCg=1, nDg=1, nDgStress=3, nCoords=2, checkZeros=True):
        """
        Initialise basic internal variables, but we must call make_cartesian_grid or make_geographic_grid before continuing.

        :param fname: Name of the file to write the output into
        """

        # Set the file name
        self.__fname = fname

        # Set basic coordinate sizes
        self.__nCg = nCg
        self.__nDg = nDg
        self.__nDgStress = nDgStress
        self.__nCoords = nCoords

        # Do we check for zeros?
        self.__checkZeros = checkZeros

    def __init_vars__(self, nFirst, nSecond):

        # Set all arrays to the correct size
        self.__nFirst = nFirst
        self.__nSecond = nSecond

        self.mask = np.zeros((nFirst, nSecond))
        self.cice = np.zeros((nFirst, nSecond))
        self.hice = np.zeros((nFirst, nSecond))
        self.hsnow = np.zeros((nFirst, nSecond))
        self.damage = np.zeros((self.__nFirst, nSecond))
        self.uice = np.zeros((nFirst, nSecond))
        self.vice = np.zeros((nFirst, nSecond))
        self.azimuth = np.zeros((nFirst, nSecond))
        self.sss = np.zeros((nFirst, nSecond))
        self.sst = np.zeros((nFirst, nSecond))

    def __init_file(self):
        """
        Open the netCDF file for writing and create the dimensions and basic structure.
        """

        self.__ncFile = netCDF4.Dataset(self.__fname, "w", format="NETCDF4")

        structure_name = "parametric_rectangular"
        self.__ncFile.structure_name = structure_name

        self.__ncFile.createDimension("ydim", self.__nFirst)
        self.__ncFile.createDimension("xdim", self.__nSecond)
        self.__ncFile.createDimension("yvertex", self.__nFirst + 1)
        self.__ncFile.createDimension("xvertex", self.__nSecond + 1)
        self.__ncFile.createDimension("y_cg", self.__nFirst * self.__nCg + 1)
        self.__ncFile.createDimension("x_cg", self.__nSecond * self.__nCg + 1)
        self.__ncFile.createDimension("dg_comp", self.__nDg)
        self.__ncFile.createDimension("dgstress_comp", self.__nDgStress)
        self.__ncFile.createDimension("ncoords", self.__nCoords)

        self.__field_dims = ("ydim", "xdim")
        self.__coord_dims = ("yvertex", "xvertex", "ncoords")

    def make_cartesian_grid(self, nFirst, nSecond, res):
        """
        Create a cartesian grid with regular spacing at a given resolution.

        :param nFirst: Number of rows (first dimension)
        :param nSecond: Number of columns (seond dimension)
        :param res: Model resolution [m]
        """

        # Set the array dimensions and resolution
        self.__res = res

        # Initialise the arrays
        self.__init_vars__(nFirst, nSecond)

        self.__init_file()

        # Array coordinates
        x = np.zeros((self.__nFirst + 1, self.__nSecond + 1))
        y = np.zeros((self.__nFirst + 1, self.__nSecond + 1))
        for j in range(self.__nFirst + 1):
            for i in range(self.__nSecond + 1):
                x[j, i] = i * self.__res
                y[j, i] = j * self.__res

        coords = self.__ncFile.createVariable("coords", "f8", self.__coord_dims)
        coords[:, :, 0] = x
        coords[:, :, 1] = y

        px = np.zeros((self.__nFirst, self.__nSecond))
        py = np.zeros((self.__nFirst, self.__nSecond))
        for j in range(self.__nFirst):
            for i in range(self.__nSecond):
                px[j, i] = (j + 0.5) * self.__res
                py[j, i] = (i + 0.5) * self.__res

        elem_x = self.__ncFile.createVariable("x", "f8", self.__field_dims)
        elem_x[:, :] = px
        elem_y = self.__ncFile.createVariable("y", "f8", self.__field_dims)
        elem_y[:, :] = py

        grid_azimuth = self.__ncFile.createVariable("grid_azimuth", "f8", self.__field_dims)
        grid_azimuth[:, :] = 0.

    def __testFields__(self):
        """
        Check if arrays are non-zero and the right size. Print a warning if
        they're zero (this may be ok). Raise a RuntimeError if the shape is wrong.
        """
        if (self.mask == 0).all():
            print("Error: 'mask' is not set (all values are zero, meaning land everywhere)")
            raise RuntimeError("'mask' is not set")

        for check in [["cice", (self.cice == 0).all(), self.cice.shape == (self.__nFirst, self.__nSecond)],
                      ["hice", (self.hice == 0).all(), self.hice.shape == (self.__nFirst, self.__nSecond)],
                      ["hsnow", (self.hsnow == 0).all(), self.hsnow.shape == (self.__nFirst, self.__nSecond)],
                      ["damage", (self.damage == 0).all(), self.damage.shape == (self.__nFirst, self.__nSecond)],
                      ["uice", (self.uice == 0).all(), self.uice.shape == (self.__nFirst, self.__nSecond)],
                      ["vice", (self.vice == 0).all(), self.vice.shape == (self.__nFirst, self.__nSecond)],
                      ["sss", (self.sss == 0).all(), self.sss.shape == (self.__nFirst, self.__nSecond)],
                      ["sst", (self.sst == 0).all(), self.sst.shape == (self.__nFirst, self.__nSecond)],
                      ["azimuth", (self.azimuth == 0).all(), self.azimuth.shape == (self.__nFirst, self.__nSecond)]]:

            if self.__checkZeros and check[1]:
                print("Warning: '" + check[0] + "' is all zeros (this may be ok, if that's what you want).")

            if not check[2]:
                print("Error: '" + check[0] + "' is the wrong shape")
                raise RuntimeError("Incorrect array shape")

    def __del__(self):
        self.__writeFile__()

    def __writeFile__(self):
        """
        Write everything to file. This is called by the destructor.
        """

        print("Producing file", self.__fname)

        self.__testFields__()

        # Set the mask
        mask = self.__ncFile.createVariable("mask", "f8", self.__field_dims)
        mask[:, :] = self.mask
        antimask = 1 - mask[:, :]

        # Set the concentration
        cice = self.__ncFile.createVariable("cice", "f8", self.__field_dims)
        cice[:, :] = self.cice

        # Set the thickness
        hice = self.__ncFile.createVariable("hice", "f8", self.__field_dims)
        hice[:, :] = self.hice

        # Set snow thickness
        hsnow = self.__ncFile.createVariable("hsnow", "f8", self.__field_dims)
        hsnow[:, :] = self.hsnow

        # Set snow thickness
        damage = self.__ncFile.createVariable("damage", "f8", self.__field_dims)
        damage[:, :] = self.damage

        # Set ice velocity
        u = self.__ncFile.createVariable("u", "f8", self.__field_dims)
        u[:, :] = self.uice
        v = self.__ncFile.createVariable("v", "f8", self.__field_dims)
        v[:, :] = self.vice

        # Set ocean state
        sst = self.__ncFile.createVariable("sst", "f8", self.__field_dims)
        sst[:, :] = self.sst
        sss = self.__ncFile.createVariable("sss", "f8", self.__field_dims)
        sss[:, :] = self.sss

        # mask data
        # TODO: Figure out why the missing data flag doesn't work
        # We should be able to define a missing value like this
        # mdi = -3.282346e38  # Minus float max
        # and then mask like this
        # cice[:, :] = cice[:, :] * mask[:, :] + antimask * mdi
        # cice.missing_value = mdi
        # but the model doesn't like this
        cice[:, :] = cice[:, :] * mask[:, :]
        hice[:, :] = hice[:, :] * mask[:, :]
        hsnow[:, :] = hsnow[:, :] * mask[:, :]
        damage[:, :] = damage[:, :] * mask[:, :]
        u[:, :] = u[:, :] * mask[:, :]
        v[:, :] = v[:, :] * mask[:, :]
        sss[:, :] = sss[:, :] * mask[:, :]
        sst[:, :] = sst[:, :] * mask[:, :]

        self.__ncFile.close()
