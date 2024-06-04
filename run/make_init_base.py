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
     1. Create an initialiser object with given filename, dimensions, and
        resolution, e.g.
      >>> init = initMaker("test", 128, 128, 1, 3e3)
     2. Modify any variables needed, e.g.
      >>> init.cice = 1
      >>> init.hice = 3
      >>> init.hice[10:20, 10:20] = 4
     3. The init file is written when the initialiser object goes out of scope
        (e.g. when the program ends or at the end of a loop).
    """

    def __init__(self, fname, nFirst, nSecond, nLayers, res, nCg=1, nDg=1, nDgStress=3, nCoords=2):
        """
        Initialise all internal variables, except __nfirst, __nsecond, __nLayers,
        and __res to zero. All arrays are set to the right size as well.
        
        :param fname: Name of the file to write the output into
        :param nFirst: Number of rows (first dimension)
        :param nSecond: Number of columns (seond dimension)
        :param nLayers: Number of thermodynamics layers (third dimension)
        :param res: Model resolution [m]
        """

        # Set the file name
        self.__fname = fname

        # Set the array dimensions and resolution
        self.__nFirst = nFirst
        self.__nSecond = nSecond
        self.__nLayers = nLayers
        self.__res = res

        # Set all arrays to the correct size
        self.mask = np.zeros((self.__nFirst, self.__nSecond))
        self.cice = np.zeros((self.__nFirst, self.__nSecond))
        self.hice = np.zeros((self.__nFirst, self.__nSecond))
        self.hsnow = np.zeros((self.__nFirst, self.__nSecond))
        self.uice = np.zeros((self.__nFirst, self.__nSecond))
        self.vice = np.zeros((self.__nFirst, self.__nSecond))
        self.azimuth = np.zeros((self.__nFirst, self.__nSecond))
        self.sss = np.zeros((self.__nFirst, self.__nSecond))
        self.sst = np.zeros((self.__nFirst, self.__nSecond))
        self.tice = np.zeros((self.__nLayers, self.__nFirst, self.__nSecond))

        # Set basic coordinate sizes
        self.__nCg = nCg
        self.__nDg = nDg
        self.__nDgStress = nDgStress
        self.__nCoords = nCoords

    def __testFields__(self):
        """
        Check if arrays are non-zero and the right size. Print a warning if
        they're zero (this may be ok). Raise a RuntimeError if the shape is wrong.
        """
        if (self.mask==0).all():
            print("Error: 'mask' is not set (all values are zero, meaning land everywhere)")
            raise RuntimeError("'mask' is not set")

        for check in [["cice", (self.cice==0).all(), self.cice.shape==(self.__nFirst,self.__nSecond)],
                      ["hice", (self.hice==0).all(), self.hice.shape==(self.__nFirst,self.__nSecond)],
                      ["hsnow", (self.hsnow==0).all(), self.hsnow.shape==(self.__nFirst,self.__nSecond)],
                      ["tice", (self.tice==0).all(), self.tice.shape==(self.__nLayers,self.__nFirst,self.__nSecond)],
                      ["uice", (self.uice==0).all(), self.uice.shape==(self.__nFirst,self.__nSecond)],
                      ["vice", (self.vice==0).all(), self.vice.shape==(self.__nFirst,self.__nSecond)],
                      ["sss", (self.sss==0).all(), self.sss.shape==(self.__nFirst,self.__nSecond)],
                      ["sst", (self.sst==0).all(), self.sst.shape==(self.__nFirst,self.__nSecond)],
                      ["azimuth", (self.azimuth==0).all(), self.azimuth.shape==(self.__nFirst,self.__nSecond)]]:

            if check[1]:
                print("Warning: '"+check[0]+"' is all zeros (this may be ok, if that's what you want).")

            if not check[2]:
                print("Error: '"+check[0]+"' is the wrong shape")
                raise RuntimeError("Incorrect array shape")

    def __del__(self):
        self.__writeFile__()

    def __writeFile__(self):
        """
        Write everything to file. This is called by the destructor.
        """

        print("Producing file", self.__fname)

        self.__testFields__()

        root = netCDF4.Dataset(self.__fname, "w", format="NETCDF4")

        structure_name = "parametric_rectangular"
        structgrp = root.createGroup("structure")
        structgrp.type = structure_name

        metagrp = root.createGroup("metadata")
        metagrp.type = structure_name
        metagrp.createGroup("configuration")  # But add nothing to it
        datagrp = root.createGroup("data")

        datagrp.createDimension("zdim", self.__nLayers)
        datagrp.createDimension("ydim", self.__nFirst)
        datagrp.createDimension("xdim", self.__nSecond)
        datagrp.createDimension("yvertex", self.__nFirst + 1)
        datagrp.createDimension("xvertex", self.__nSecond + 1)
        datagrp.createDimension("y_cg", self.__nFirst * self.__nCg + 1)
        datagrp.createDimension("x_cg", self.__nSecond * self.__nCg + 1)
        datagrp.createDimension("dg_comp", self.__nDg)
        datagrp.createDimension("dgstress_comp", self.__nDgStress)
        datagrp.createDimension("ncoords", self.__nCoords)

        field_dims = ("ydim", "xdim")
        coord_dims = ("yvertex", "xvertex", "ncoords")

        # Array coordinates
        x = np.zeros((self.__nFirst + 1, self.__nSecond + 1))
        y = np.zeros((self.__nFirst + 1, self.__nSecond + 1))
        for j in range(self.__nFirst + 1):
            for i in range(self.__nSecond + 1):
                x[j, i] = i * self.__res
                y[j, i] = j * self.__res

        coords = datagrp.createVariable("coords", "f8", coord_dims)
        coords[:, :, 0] = x
        coords[:, :, 1] = y

        px = np.zeros((self.__nFirst, self.__nSecond))
        py = np.zeros((self.__nFirst, self.__nSecond))
        for j in range(self.__nFirst):
            for i in range(self.__nSecond):
                px[j, i] = (j + 0.5) * self.__res
                py[j, i] = (i + 0.5) * self.__res

        elem_x = datagrp.createVariable("x", "f8", field_dims)
        elem_x[:, :] = px
        elem_y = datagrp.createVariable("y", "f8", field_dims)
        elem_y[:, :] = py

        grid_azimuth = datagrp.createVariable("grid_azimuth", "f8", field_dims)
        grid_azimuth[:, :] = self.azimuth

        # Set the mask
        mask = datagrp.createVariable("mask", "f8", field_dims)
        mask[:, :] = self.mask
        antimask = 1 - mask[:, :]

        # Set the concentration
        cice = datagrp.createVariable("cice", "f8", field_dims)
        cice[:, :] = self.cice

        # Set the thickness
        hice = datagrp.createVariable("hice", "f8", field_dims)
        hice[:, :] = self.hice

        # Set snow thickness
        hsnow = datagrp.createVariable("hsnow", "f8", field_dims)
        hsnow[:, :] = self.hsnow

        # Set ice temperatures
        tice = datagrp.createVariable("tice", "f8", ("zdim", "ydim", "xdim"))
        tice[:, :, :] = self.tice

        # Set ice velocity
        u = datagrp.createVariable("u", "f8", field_dims)
        u[:, :] = self.uice
        v = datagrp.createVariable("v", "f8", field_dims)
        v[:, :] = self.vice

        # Set ocean state
        sst = datagrp.createVariable("sst", "f8", field_dims)
        sst[:, :] = self.sst
        sss = datagrp.createVariable("sss", "f8", field_dims)
        sss[:, :] = self.sss

        # mask data
        mdi = -3.282346e38  # Minus float max
        cice[:, :] = cice[:, :] * mask[:, :] + antimask * mdi
        cice.missing_value = mdi
        hice[:, :] = hice[:, :] * mask[:, :] + antimask * mdi
        hice.missing_value = mdi
        hsnow[:, :] = hsnow[:, :] * mask[:, :] + antimask * mdi
        hsnow.missing_value = mdi
        u[:, :] = u[:, :] * mask[:, :] + antimask * mdi
        u.missing_value = mdi
        v[:, :] = v[:, :] * mask[:, :] + antimask * mdi
        v.missing_value = mdi
        tice[:, :, :] = tice[:, :, :] * mask[:, :] + antimask * mdi
        tice.missing_value = mdi
        grid_azimuth[:, :] = grid_azimuth[:, :] * mask[:, :] + antimask * mdi
        grid_azimuth.missing_value = mdi
        sss[:, :] = sss[:, :] * mask[:, :] + antimask * mdi
        sss.missing_value = mdi
        sst[:, :] = sst[:, :] * mask[:, :] + antimask * mdi
        sst.missing_value = mdi

        root.close()
