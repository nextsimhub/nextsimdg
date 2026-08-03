import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import Delaunay


class FloodInterpolator:
    """
    Flood land cells using a precomputed Delaunay interpolation operator.

    The expensive geometric preprocessing is performed once during
    construction and can then be reused to flood any number of fields.
    """

    def __init__(self, x, y, land_mask):
        """
        Initialise an object containing the Delaunay triangulation used for flooding.

        Example usage:

        >> from flooding import FloodInterpolator
        >> flooder = FloodInterpolator(ds.x.values, ds.y.values, ds.mask.values.astype(bool))

        where ds is a xarray Dataset object containing the selected variables.

        Parameters
        ----------
        x : (nx,) array
            x coordinates
        y : (ny,) array
            y coordinates
        land_mask : (ny,nx) bool
            True over land, False over ocean.
        """
        self.A, self.ocean, self.land, self.valid = self.__build_flood_operator(
            x, y, land_mask
        )
        self.land_mask = land_mask

    def __build_flood_operator(self, x, y, land_mask):
        """
        Construct a sparse interpolation operator that floods land cells from neighbouring ocean cells using linear interpolation.

        Parameters
        ----------
        x : (nx,) array
            x coordinates
        y : (ny,) array
            y coordinates
        land_mask : (ny,nx) bool
            True over land, False over ocean.

        Returns
        -------
        A : csr_matrix
            Sparse interpolation matrix.
        ocean : bool array
            Ocean mask.
        land : bool array
            Land mask.
        valid : bool array
            Land points inside the convex hull.
        """
        ocean = ~land_mask
        land = land_mask

        X, Y = np.meshgrid(x, y)

        pts = np.c_[X[ocean], Y[ocean]]
        targets = np.c_[X[land], Y[land]]

        tri = Delaunay(pts)

        simplex = tri.find_simplex(targets)
        valid = simplex >= 0

        vertices = tri.simplices[simplex[valid]]

        T = tri.transform[simplex[valid]]
        delta = targets[valid] - T[:, 2]

        bary = np.einsum("ijk,ik->ij", T[:, :2], delta)

        weights = np.column_stack((bary, 1.0 - bary.sum(axis=1)))

        rows = np.repeat(np.arange(valid.sum()), 3)

        A = csr_matrix(
            (weights.ravel(), (rows, vertices.ravel())),
            shape=(valid.sum(), ocean.sum()),
        )

        return A, ocean, land, valid

    def flood_field(self, data):
        """
        Flood a field using a precomputed sparse interpolation operator.

        Parameters
        ----------
        data : ndarray
            Array with shape (..., ny, nx). The last two dimensions must
            correspond to the grid used to build the interpolation operator.
            The input array is modified in place.
        """
        nland = self.land.sum()

        # Choose an appropriate floating-point dtype
        dtype = np.result_type(data.dtype, np.float32)
        filled_land = np.empty(nland, dtype=dtype)

        if data.ndim == 2:
            filled_land.fill(np.nan)
            filled_land[self.valid] = self.A @ data[self.ocean]
            data[self.land] = filled_land

        else:
            for index in np.ndindex(data.shape[:-2]):
                field = data[index]

                filled_land.fill(np.nan)
                filled_land[self.valid] = self.A @ field[self.ocean]

                field[self.land] = filled_land

        return data

    def flood_dataset(self, ds, fields="all"):
        """
        Flood selected variables of a xarray dataset.

        Parameters
        ----------
        ds : xarray dataset
            The dataset to be flooded.
        fields : iterable of str
            A list of the variable names to be flooded, or the keyword "all" to
            process all fields which have more than two dimensions.
        """
        my_fields = []
        if fields == "all":
            for var in list(ds.variables):
                if len(ds.variables[var].dims) > 2:
                    my_fields.append(var)
        else:
            my_fields = fields

        for var in my_fields:
            self.flood_field(ds[var].values)

    def flood_zeros(self, ds, fields):
        """
        Flood selected variables of a xarray dataset with zeros.

        Parameters
        ----------
        ds : xarray dataset
            The dataset to be flooded.
        fields : iterable of str
            A list of the variable names to be flooded.
        """
        for var in fields:
            ds[var].values[np.stack([self.land_mask] * ds[var].values.shape[0])] = 0
