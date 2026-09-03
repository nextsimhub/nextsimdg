"""
Area-average IBCAO bathymetry onto a rotated rectangular model grid.

The IBCAO GeoTIFF and the model grid are both defined in EPSG:3996.
The model grid may be arbitrarily rotated relative to the native IBCAO
grid.

IBCAO convention used here:

    depth < 0    -> ocean
    depth >= 0   -> land

The main operation is:

    IBCAO 400 m raster
            |
            | area-weighted averaging
            v
    rotated model grid

average_to_grid() returns:

    depth
        Mean ocean depth within each model cell [m].

    ocean_mask
        Boolean mask indicating which model cells are ocean [True, False].
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from projection import PolarStereoGrid, PolarStereoProjection
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.warp import reproject


class IBCAO:
    """
    IBCAO bathymetry reader and area-averager.

    Parameters
    ----------
    filename : str or pathlib.Path
        IBCAO GeoTIFF.
    """

    def __init__(self, filename):

        self.filename = Path(filename)

        if not self.filename.exists():
            msg = f"IBCAO file not found: {self.filename}"
            raise FileNotFoundError(msg)

        self.src = rasterio.open(self.filename)

        # ------------------------------------------------------------------
        # CRS
        # ------------------------------------------------------------------

        self.crs = self.src.crs

        if self.crs is None:
            msg = "IBCAO GeoTIFF does not contain a CRS"
            raise ValueError(msg)

        if self.crs.to_epsg() != 3996:
            msg = f"Expected IBCAO CRS EPSG:3996, found {self.crs}"
            raise ValueError(msg)

        # ------------------------------------------------------------------
        # Raster metadata
        # ------------------------------------------------------------------

        self.transform = self.src.transform

        self.width = self.src.width
        self.height = self.src.height

        self.bounds = self.src.bounds

        self.nodata = self.src.nodata

        self.dx = abs(self.transform.a)
        self.dy = abs(self.transform.e)

    # ------------------------------------------------------------------
    # Context manager
    # ------------------------------------------------------------------

    def close(self):
        """Close the GeoTIFF."""
        self.src.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    # ------------------------------------------------------------------
    # Target transform
    # ------------------------------------------------------------------

    @staticmethod
    def _target_transform(
        grid: PolarStereoGrid, projection: PolarStereoProjection
    ) -> Affine:

        # Lower-left model-grid corner -> EPSG:3996.
        x0, y0 = projection.unrotate(grid.xbnds[0], grid.ybnds[0])

        # Model +x basis vector in EPSG:3996.
        ex_x, ex_y = projection.unrotate(grid.dx, 0.0)

        # Model +y basis vector in EPSG:3996.
        ey_x, ey_y = projection.unrotate(0.0, grid.dy)

        # Upper-left corner.
        x_ul = x0 + grid.ny * ey_x
        y_ul = y0 + grid.ny * ey_y

        return Affine(ex_x, -ey_x, x_ul, ex_y, -ey_y, y_ul)

    # ------------------------------------------------------------------
    # Area averaging
    # ------------------------------------------------------------------

    def average_to_grid(self, grid: PolarStereoGrid, projection: PolarStereoProjection):
        """
        Area-average IBCAO bathymetry onto the model grid.

        Parameters
        ----------
        grid : PolarStereoGrid
            Target model grid.
        projection : PolarStereoProjection
            Projection defining the model-grid rotation.

        Returns
        -------
        depth : ndarray
            Mean ocean depth in each model cell [m].

            NaN is returned for cells containing no ocean.

        ocean_mask : ndarray
            Boolean mask indicating which model cells are ocean [True, False].
        """
        target_shape = (grid.ny, grid.nx)
        target_transform = self._target_transform(grid, projection)

        # ------------------------------------------------------------------
        # Read IBCAO
        # ------------------------------------------------------------------
        source = self.src.read(1).astype(np.float32, copy=False)

        # ------------------------------------------------------------------
        # Elevation field
        # ------------------------------------------------------------------
        # We don't set land to zero, so that the interpolation respects the coast line
        ocean_depth = source.astype(np.float32)

        # ------------------------------------------------------------------
        # Resample depth
        # ------------------------------------------------------------------
        depth_sum = np.zeros(target_shape, dtype=np.float32)
        reproject(
            source=ocean_depth,
            destination=depth_sum,
            src_transform=self.transform,
            src_crs=self.crs,
            dst_transform=target_transform,
            dst_crs=self.crs,
            resampling=Resampling.average,
        )
        depth_sum = np.flipud(depth_sum)

        # ------------------------------------------------------------------
        # Set depth to NaN for cells with no ocean and define an ocean mask
        # ------------------------------------------------------------------
        depth = np.full(target_shape, np.nan, dtype=np.float32)

        ocean_cells = depth_sum < 0.0
        depth[ocean_cells] = depth_sum[ocean_cells]

        return depth, ocean_cells
