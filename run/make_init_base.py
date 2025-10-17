import netCDF4
import numpy as np


class initMaker:
    """
    A "plug-and-play" initialisation class for neXtSIM.

    The user needs to supply at minimum the grid dimensions and resolution. They may
    also supply any initialisation fields they need, as well as a land mask.

    Usage:
     0. Import make_init_base
      >> from make_init_base import initMaker
     1. Create an initialiser object with a given filename, e.g.
      >> init = initMaker("test")
         2a. Create a cartesian grid from dimensions and resolution, e.g.
          >> init.make_cartesian_grid(128, 128, 3e3)
         2b. Create a geographic grid from a file containing coordinates. We need to define type of grid:
            * "p": only centre points of the grid are given
            * "ur": the centre and upper right corner of the grid is given
            * "ll": the centre and lower left corner of the grid is given
          init.make_geographic_grid("coordinates.nc", "ur", plat_name="gphit", plon_name="glamt", qlat_name="gphif", qlon_name="glamf")
     3. Modify any variables needed, e.g.
      >> init.cice = 1
      >> init.hice = 3
      >> init.hice[10:20, 10:20] = 4
     4. The init file is written when the initialiser object goes out of scope (e.g. when the program ends or at the end of a loop).
    """

    def __init__(
        self,
        fname,
        nCg=1,
        nDg=1,
        nDgStress=3,
        nCoords=2,
        checkZeros=True,
    ):
        """
        Initialise basic internal variables, but we must call make_cartesian_grid or make_geographic_grid before continuing.

        :param fname: Name of the file to write the output into
        """
        # Set the file name
        self._fname = fname

        # Set basic coordinate sizes
        self._nCg = nCg
        self._nDg = nDg
        self._nDgStress = nDgStress
        self._nCoords = nCoords

        # Do we check for zeros?
        self._checkZeros = checkZeros

    def _init_vars_and_file(self, nFirst, nSecond):
        """
        Initialise the arrays. Open the netCDF file for writing and create the dimensions and basic structure. This must
        be called from inside make_cartesian_grid or make_geographic_grid.

        :param nFirst: The number of rows (first dimension)
        :param nSecond: The number of columns (second dimension)
        """
        # Set all arrays to the correct size
        self._nFirst = nFirst
        self._nSecond = nSecond

        # Set all arrays to zero ...
        self.mask = np.zeros((nFirst, nSecond))
        self.cice = np.zeros((nFirst, nSecond))
        self.hice = np.zeros((nFirst, nSecond))
        self.hsnow = np.zeros((nFirst, nSecond))
        self.uice = np.zeros((nFirst, nSecond))
        self.vice = np.zeros((nFirst, nSecond))
        self.azimuth = np.zeros((nFirst, nSecond))
        self.sss = np.zeros((nFirst, nSecond))
        self.sst = np.zeros((nFirst, nSecond))

        # ... except for the damage array, which is all ones.
        self.damage = np.ones((nFirst, nSecond))

        # Open the netCDF file for writing and create the dimensions and basic structure.
        self._ncFile = netCDF4.Dataset(self._fname, "w", format="NETCDF4")

        structure_name = "parametric_rectangular"
        self._ncFile.structure_name = structure_name

        self._ncFile.createDimension("ydim", nFirst)
        self._ncFile.createDimension("xdim", nSecond)
        self._ncFile.createDimension("yvertex", nFirst + 1)
        self._ncFile.createDimension("xvertex", nSecond + 1)
        self._ncFile.createDimension("y_cg", nFirst * self._nCg + 1)
        self._ncFile.createDimension("x_cg", nSecond * self._nCg + 1)
        self._ncFile.createDimension("dg_comp", self._nDg)
        self._ncFile.createDimension("dgstress_comp", self._nDgStress)
        self._ncFile.createDimension("ncoords", self._nCoords)

        self._field_dims = ("ydim", "xdim")
        self._coord_dims = ("yvertex", "xvertex", "ncoords")

    def make_cartesian_grid(self, nFirst, nSecond, res):
        """
        Create a cartesian grid with regular spacing at a given resolution.

        :param nFirst: Number of rows (first dimension)
        :param nSecond: Number of columns (second dimension)
        :param res: Model resolution [m]
        """
        # Initialise the arrays and the output file
        self._init_vars_and_file(nFirst, nSecond)

        # Array coordinates
        x = np.zeros((nFirst + 1, nSecond + 1))
        y = np.zeros((nFirst + 1, nSecond + 1))
        for j in range(nFirst + 1):
            for i in range(nSecond + 1):
                x[j, i] = i * res
                y[j, i] = j * res

        coords = self._ncFile.createVariable("coords", "f8", self._coord_dims)
        coords[:, :, 0] = x
        coords[:, :, 1] = y

        px = np.zeros((nFirst, nSecond))
        py = np.zeros((nFirst, nSecond))
        for j in range(nFirst):
            for i in range(nSecond):
                px[j, i] = (j + 0.5) * res
                py[j, i] = (i + 0.5) * res

        elem_x = self._ncFile.createVariable("x", "f8", self._field_dims)
        elem_x[:, :] = px
        elem_y = self._ncFile.createVariable("y", "f8", self._field_dims)
        elem_y[:, :] = py

        grid_azimuth = self._ncFile.createVariable(
            "grid_azimuth", "f8", self._field_dims
        )
        grid_azimuth[:, :] = 0.0

    def make_geographic_grid(
        self,
        grid_file,
        pos,
        plon_name="plon",
        plat_name="plat",
        qlon_name="qlon",
        qlat_name="qlat",
    ):
        """
        Create a geographic grid from a file containing coordinates.

        :param grid_file: The input grid file.
        :param pos: The position of the grid corners. Either "ur", "ll", or "p".
        :param plon_name: Name of the longitude variable in the grid file for the grid cell centre.
        :param plat_name: Name of the latitude variable in the grid file for the grid cell centre.
        :param qlon_name: Name of the longitude variable in the grid file for the grid cell upper right or lower left corner.
        :param qlat_name: Name of the latitude variable in the grid file for the grid cell upper right or lower left corner.
        :return:
        """
        grid = netCDF4.Dataset(f"{grid_file}", "r")

        # Grid dimensions. We're dealing with files written in and for FORTRAN, so y is the first dimension and x the second.
        nfirst = grid.dimensions["y"].size
        nsecond = grid.dimensions["x"].size

        # Initialise the arrays and the output file
        self._init_vars_and_file(nfirst, nsecond)

        # Array coordinates
        node_lon = np.zeros((nfirst + 1, nsecond + 1))
        node_lat = np.zeros((nfirst + 1, nsecond + 1))

        # Deduce the coordinates of the grid corners
        if pos == "ur":
            """
            This is for grids with plon and plat at the centre of the grid cell and qlon and qlat at the upper right
            corner; typical for an ocean model C-grid. The NEMO/ORCA grid is one of those.
            That means making up the locations of the first few points.
            """
            self._plon = grid.variables[plon_name][:, :]
            self._plat = grid.variables[plat_name][:, :]
            qlon = grid.variables[qlon_name][:, :]
            qlat = grid.variables[qlat_name][:, :]

            # Stay in the [-180, 180] range
            self._plon = self._wrap_to_180(self._plon)

            # We proceed through the grid from the lower left corner in a row-major order. For the UR grid, that means
            # making up the locations of the first few points.

            # Extrapolate to get the lower left corner
            node_lon[0, 0] = (
                qlon[0, 0] - (qlon[1, 0] - qlon[0, 0]) - (qlon[0, 1] - qlon[0, 0])
            )
            node_lat[0, 0] = (
                qlat[0, 0] - (qlat[1, 0] - qlat[0, 0]) - (qlat[0, 1] - qlat[0, 0])
            )

            # Extrapolate to get the rest of the bottom row
            node_lon[0, 1:] = qlon[0, :] - (qlon[1, :] - qlon[0, :])
            node_lat[0, 1:] = qlat[0, :] - (qlat[1, :] - qlat[0, :])

            # Extrapolate to get the first column
            node_lon[1:, 0] = qlon[:, 0] - (qlon[:, 1] - qlon[:, 0])
            node_lat[1:, 0] = qlat[:, 0] - (qlat[:, 1] - qlat[:, 0])

            # Fill the rest with the values provided on the upper right corner
            node_lon[1:, 1:] = qlon[:, :]
            node_lat[1:, 1:] = qlat[:, :]

        elif pos == "p":
            """
            This is for grids with plon and plat at the centre of the grid cell.
            That means making up the locations of the grid nodes/vertices.
            """
            self._plon = grid.variables[plon_name][:, :]
            self._plat = grid.variables[plat_name][:, :]

            # Stay in the [-180, 180] range
            self._plon = self._wrap_to_180(self._plon)

            # First, we project the grid centers to Greenland
            (rlat, rlon) = self._rotate_pole_to_greenland(self._plat, self._plon)

            # Then, we create the grid nodes/vertices

            # Differences in longitude and latitude
            dlat_dx = 0.5 * np.diff(rlat, axis=0)
            dlon_dx = 0.5 * self._wrap_to_180(np.diff(rlon, axis=0))

            # A minus here, because the array is effectively flipped upside down.
            dlat_dy = -0.5 * np.diff(rlat, axis=1)
            dlon_dy = -0.5 * self._wrap_to_180(np.diff(rlon, axis=1))

            # The corners
            # lower left
            node_lon[0, 0] = rlon[0, 0] - dlon_dx[0, 0] - dlon_dy[0, 0]
            node_lat[0, 0] = rlat[0, 0] - dlat_dx[0, 0] - dlat_dy[0, 0]
            # upper left
            node_lon[0, -1] = rlon[0, -1] - dlon_dx[0, -1] + dlon_dy[0, -1]
            node_lat[0, -1] = rlat[0, -1] - dlat_dx[0, -1] + dlat_dy[0, -1]

            # upper right
            node_lon[-1, -1] = rlon[-1, -1] + dlon_dx[-1, -1] + dlon_dy[-1, -1]
            node_lat[-1, -1] = rlat[-1, -1] + dlat_dx[-1, -1] + dlat_dy[-1, -1]
            # lower right
            node_lon[-1, 0] = rlon[-1, 0] + dlon_dx[-1, 0] - dlon_dy[-1, 0]
            node_lat[-1, 0] = rlat[-1, 0] + dlat_dx[-1, 0] - dlat_dy[-1, 0]

            # Columns and rows
            # first column
            node_lon[0, 1:-1] = rlon[0, :-1] - dlon_dx[0, :-1] + dlon_dy[0, :]
            node_lat[0, 1:-1] = rlat[0, :-1] - dlat_dx[0, :-1] + dlat_dy[0, :]
            # top row
            node_lon[1:-1, -1] = rlon[:-1, -1] + dlon_dx[:, -1] + dlon_dy[:-1, -1]
            node_lat[1:-1, -1] = rlat[:-1, -1] + dlat_dx[:, -1] + dlat_dy[:-1, -1]
            # last column
            node_lon[-1, 1:-1] = rlon[-1, :-1] + dlon_dx[-1, :-1] + dlon_dy[-1, :]
            node_lat[-1, 1:-1] = rlat[-1, :-1] + dlat_dx[-1, :-1] + dlat_dy[-1, :]
            # bottom row
            node_lon[1:-1, 0] = rlon[:-1, 0] + dlon_dx[:, 0] - dlon_dy[:-1, 0]
            node_lat[1:-1, 0] = rlat[:-1, 0] + dlat_dx[:, 0] - dlat_dy[:-1, 0]

            # Fill the innermost points
            # The upper right corner is centre coordinate, plus half the grid spacing in each direction
            node_lon[1:-1, 1:-1] = rlon[:-1, :-1] + dlon_dx[:, 1:] + dlon_dy[1:, :]
            node_lat[1:-1, 1:-1] = rlat[:-1, :-1] + dlat_dx[:, 1:] + dlat_dy[1:, :]

            (node_lat, node_lon) = self._rotate_pole_from_greenland(node_lat, node_lon)
            node_lon = self._wrap_to_180(node_lon)

        else:
            raise ValueError(
                f"Position {pos} not yet implemented (expected 'ur', 'll', or 'p')."
            )

        coords = self._ncFile.createVariable("coords", "f8", self._coord_dims)
        coords[:, :, 0] = node_lon
        coords[:, :, 1] = node_lat

        elem_lon = self._ncFile.createVariable("longitude", "f8", self._field_dims)
        elem_lon[:, :] = self._plon
        elem_lat = self._ncFile.createVariable("latitude", "f8", self._field_dims)
        elem_lat[:, :] = self._plat

        grid_azimuth = self._ncFile.createVariable(
            "grid_azimuth", "f8", self._field_dims
        )

        # Rotate the pole to Greenland
        (rlat, rlon) = self._rotate_pole_to_greenland(node_lat, node_lon)
        # Return the grid azimuth
        grid_azimuth[:, :] = self._grid_angle(rlon, rlat)

        grid.close()

    def _rotate_pole_to_greenland(self, lat, lon):
        """
        Rotates the mesh such that the singularities are in Greenland / Antarctica at 75°N / 40°W and 75°S / 140°E
        This is a copy of ParametricMesh::RotatePoleToGreenland in dynamics/src/include/ParametricMesh.hpp.

        :param lat: Latitudes of the mesh
        :param lon: Longitudes of the mesh
        :return: A tuple of the rotated latitudes and longitudes.
        """
        latr = np.deg2rad(lat)
        lonr = np.deg2rad(lon)

        x = np.cos(latr) * np.cos(lonr)
        y = np.cos(latr) * np.sin(lonr)
        z = np.sin(latr)

        aw = np.deg2rad(40.0)
        x1 = np.cos(aw) * x - np.sin(aw) * y
        y1 = np.sin(aw) * x + np.cos(aw) * y
        z1 = z

        bw = np.deg2rad(15.0)
        x2 = np.cos(bw) * x1 - np.sin(bw) * z1
        y2 = y1
        z2 = np.sin(bw) * x1 + np.cos(bw) * z1

        return (np.rad2deg(np.arcsin(z2)), np.rad2deg(np.arctan2(y2, x2)))

    def _rotate_pole_from_greenland(self, lat, lon):
        """
        Rotates back the mesh from having the singularities in Greenland / Antarctica at 75°N / 40°W and 75°S / 140°E to the original coordinates.
        This is a copy of ParametricMesh::RotatePoleFromGreenland in dynamics/src/include/ParametricMesh.hpp.

        :param lat: Latitudes of the rotated mesh
        :param lon: Longitudes of the rotated mesh
        :return: A tuple of the geographical latitudes and longitudes.
        """
        latr = np.deg2rad(lat)
        lonr = np.deg2rad(lon)

        x = np.cos(latr) * np.cos(lonr)
        y = np.cos(latr) * np.sin(lonr)
        z = np.sin(latr)

        aw = np.deg2rad(-40.0)
        bw = np.deg2rad(-15.0)

        x1 = np.cos(bw) * x - np.sin(bw) * z
        y1 = y
        z1 = np.sin(bw) * x + np.cos(bw) * z

        x2 = np.cos(aw) * x1 - np.sin(aw) * y1
        y2 = np.sin(aw) * x1 + np.cos(aw) * y1
        z2 = z1

        return (np.rad2deg(np.arcsin(z2)), np.rad2deg(np.arctan2(y2, x2)))

    def _wrap_to_180(self, x_in):
        """
        Wrap the input array to the range [-180, 180].

        :param x_in: Input longitudes.
        :return: Output longitudes in the range [-180, 180].
        """
        x = x_in
        x += 180.0
        x %= 360.0
        x -= 180.0
        return x

    def _grid_angle(self, lons, lats):
        """
        Calculate the angle between the grid y-axis and north.

        :param lons: Longitudes of the grid vertices.
        :param lats: Latitudes of the grid vertices.
        :return: Angle between the grid y-axis and north.
        """
        n = np.shape(lons)[0] - 1
        m = np.shape(lons)[1] - 1

        theta = np.zeros((n, m))

        for i in range(n):
            for j in range(m):
                lonVerts = lons[i : i + 2, j : j + 2]
                latVerts = lats[i : i + 2, j : j + 2]

                # Rotate the longitude if needed
                if (
                    abs(lonVerts[0, 0] - lonVerts[1, 0]) > 90.0
                    or abs(lonVerts[0, 1] - lonVerts[1, 1]) > 90
                    or abs(lonVerts[0, 0] - lonVerts[0, 1]) > 90
                    or abs(lonVerts[1, 0] - lonVerts[1, 1]) > 90
                ):
                    for ii in range(2):
                        for jj in range(2):
                            if lonVerts[ii, jj] < 0:
                                lonVerts[ii, jj] += 360

                # Find the centre points along the x-axis
                # The array looks like this:
                # [0,0] --- (lat0,lon0) --- [0,1]
                #   |                         |
                #   |                         |
                #   |                         |
                # [1,0] --- (lat1,lon1) --- [1,1]
                # Where (lat1,lon1) and (lat0,lon0) are the centre points of the top and bottom sides of the grid cell. I've
                # removed the division by two from the averaging, as it gets cancelled out in atan2.
                lat0 = latVerts[0, 0] + latVerts[0, 1]
                lat1 = latVerts[1, 0] + latVerts[1, 1]
                lon0 = lonVerts[0, 0] + lonVerts[0, 1]
                lon1 = lonVerts[1, 0] + lonVerts[1, 1]

                # Calculate the angle as the change in lat/lat coordinates across the grid cell
                dx = lat1 - lat0
                dy = lon1 - lon0

                # Negative of the angle to get the rotation of the vector to the displaced pole (not the other way around).
                # Negative on dy, because the y-axis is pointing downwards.
                theta[i, j] = -np.arctan2(-dy, dx)

        return theta

    def get_element_longitude(self):
        """
        Returns the longitude of the grid cell centre.

        :return: Longitude of the grid cell centre.
        """
        return self._plon

    def get_element_latitude(self):
        """
        Returns the latitude of the grid cell centre.

        :return: Latitude of the grid cell centre.
        """
        return self._plat

    def _testFields(self):
        """
        Check if arrays are non-zero and the right size.

        Print a warning if they're zero (this may be ok). Raise a RuntimeError if the
        shape is wrong.
        """
        if (self.mask == 0).all():
            print(
                "Error: 'mask' is not set (all values are zero, meaning land everywhere)"
            )
            runtime_err = "'mask' is not set"
            raise RuntimeError(runtime_err)

        for field, allZero, wrongShape in [
            [
                "cice",
                (self.cice == 0).all(),
                self.cice.shape == (self._nFirst, self._nSecond),
            ],
            [
                "hice",
                (self.hice == 0).all(),
                self.hice.shape == (self._nFirst, self._nSecond),
            ],
            [
                "hsnow",
                (self.hsnow == 0).all(),
                self.hsnow.shape == (self._nFirst, self._nSecond),
            ],
            [
                "sss",
                (self.sss == 0).all(),
                self.sss.shape == (self._nFirst, self._nSecond),
            ],
            [
                "sst",
                (self.sst == 0).all(),
                self.sst.shape == (self._nFirst, self._nSecond),
            ],
        ]:
            if self._checkZeros and allZero:
                print(
                    f"Warning: '{field}' is all zeros (this may be ok, if that's what you want)."
                )

            if not wrongShape:
                print(f"Error: '{field}' is the wrong shape")
                runtime_err = "Incorrect array shape"
                raise RuntimeError(runtime_err)

    def __del__(self):
        """Destructor that writes the file when the object goes out of scope."""
        self._writeFile()

    def _writeFile(self):
        """Write everything to a file. This is called by the destructor."""
        print("Producing file", self._fname)

        self._testFields()

        # Set the mask
        mask = self._ncFile.createVariable("mask", "f8", self._field_dims)
        mask[:, :] = self.mask
        antimask = 1 - mask[:, :]

        # Set the concentration
        cice = self._ncFile.createVariable("cice", "f8", self._field_dims)
        cice[:, :] = self.cice

        # Set the thickness
        hice = self._ncFile.createVariable("hice", "f8", self._field_dims)
        hice[:, :] = self.hice

        # Set snow thickness
        hsnow = self._ncFile.createVariable("hsnow", "f8", self._field_dims)
        hsnow[:, :] = self.hsnow

        # Set snow thickness
        damage = self._ncFile.createVariable("damage", "f8", self._field_dims)
        damage[:, :] = self.damage

        # Set ice velocity
        u = self._ncFile.createVariable("u", "f8", self._field_dims)
        u[:, :] = self.uice
        v = self._ncFile.createVariable("v", "f8", self._field_dims)
        v[:, :] = self.vice

        # Set ocean state
        sst = self._ncFile.createVariable("sst", "f8", self._field_dims)
        sst[:, :] = self.sst
        sss = self._ncFile.createVariable("sss", "f8", self._field_dims)
        sss[:, :] = self.sss

        # mask data
        mdi = -3.282346e38  # Minus float max

        cice[:, :] = cice[:, :] * mask[:, :] + antimask * mdi
        cice.missing_value = mdi

        hice[:, :] = hice[:, :] * mask[:, :] + antimask * mdi
        hice.missing_value = mdi

        hsnow[:, :] = hsnow[:, :] * mask[:, :] + antimask * mdi
        hsnow.missing_value = mdi

        damage[:, :] = damage[:, :] * mask[:, :] + antimask * mdi
        damage.missing_value = mdi

        # U and V must be zero on land
        u[:, :] = u[:, :] * mask[:, :]
        v[:, :] = v[:, :] * mask[:, :]

        sss[:, :] = sss[:, :] * mask[:, :] + antimask * mdi
        sss.missing_value = mdi

        sst[:, :] = sst[:, :] * mask[:, :] + antimask * mdi
        sst.missing_value = mdi

        self._ncFile.close()
