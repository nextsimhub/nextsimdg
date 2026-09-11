"""
Utilities for generating a rotated Arctic polar stereographic grid.

The base projection is EPSG:3996, the projection used by IBCAO.

The model grid is a regular rectangular grid in EPSG:3996 coordinates,
with an optional rotation applied in the projected Cartesian plane.

Coordinate convention
---------------------

Unrotated coordinates:

    EPSG:3996 x, y
          │
          │ rotate
          ▼
    model x, y

Thus, when rotation == 0, model coordinates are identical to IBCAO
coordinates.

The grid limits xmin/xmax/ymin/ymax refer to the OUTER CELL EDGES.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyproj import CRS, Transformer


@dataclass(slots=True)
class PolarStereoGrid:
    """
    Rectangular model grid.

    Parameters
    ----------
    x, y : ndarray
        1D arrays containing cell-centre coordinates [m].
    xbnds, ybnds : ndarray
        1D arrays containing cell-edge coordinates [m].
    X, Y : ndarray
        2D arrays of cell-centre coordinates [m].
    lon, lat : ndarray
        2D longitude and latitude of cell centres.
    dx, dy : float
        Grid spacing [m].
    """

    x: np.ndarray
    y: np.ndarray

    xbnds: np.ndarray
    ybnds: np.ndarray

    X: np.ndarray
    Y: np.ndarray

    lon: np.ndarray
    lat: np.ndarray

    dx: float
    dy: float

    lon_corner: np.ndarray
    lat_corner: np.ndarray

    @property
    def nx(self) -> int:
        return self.x.size

    @property
    def ny(self) -> int:
        return self.y.size

    @property
    def shape(self) -> tuple[int, int]:
        return self.X.shape

    @property
    def cell_area(self) -> float:
        """Area of one model grid cell [m²]."""
        return self.dx * self.dy


@dataclass(slots=True)
class PolarStereoProjection:
    """
    Rotated Arctic polar stereographic projection.

    The underlying CRS is EPSG:3996 (IBCAO polar stereographic).

    Parameters
    ----------
    rotation : float
        Counter-clockwise rotation angle in the projected plane [degrees].
    """

    rotation: float = 0.0
    epsg: int = 3996
    crs: CRS = None
    _to_xy: Transformer = None
    _to_lonlat: Transformer = None
    cos_theta: float = None
    sin_theta: float = None

    def __post_init__(self):

        self.crs = CRS.from_epsg(self.epsg)

        self._to_xy = Transformer.from_crs(
            CRS.from_epsg(4326), self.crs, always_xy=True
        )

        self._to_lonlat = Transformer.from_crs(
            self.crs, CRS.from_epsg(4326), always_xy=True
        )

        theta = np.deg2rad(self.rotation)

        self.cos_theta = np.cos(theta)
        self.sin_theta = np.sin(theta)

    # ------------------------------------------------------------------
    # Rotation
    # ------------------------------------------------------------------

    def rotate(self, x, y):
        """
        Rotate EPSG:3996 coordinates into model coordinates.

        The rotation is counter-clockwise.
        """
        xr = self.cos_theta * x - self.sin_theta * y
        yr = self.sin_theta * x + self.cos_theta * y

        return xr, yr

    # ------------------------------------------------------------------

    def unrotate(self, x, y):
        """Convert model coordinates back to unrotated EPSG:3996 coordinates."""
        xr = self.cos_theta * x + self.sin_theta * y
        yr = -self.sin_theta * x + self.cos_theta * y

        return xr, yr

    # ------------------------------------------------------------------
    # Geographic conversion
    # ------------------------------------------------------------------

    def lonlat_to_xy(self, lon, lat):
        """Convert longitude/latitude to rotated model coordinates."""
        x, y = self._to_xy.transform(lon, lat)

        return self.rotate(x, y)

    # ------------------------------------------------------------------

    def xy_to_lonlat(self, x, y):
        """Convert rotated model coordinates to longitude/latitude."""
        x_unrot, y_unrot = self.unrotate(x, y)

        lon, lat = self._to_lonlat.transform(x_unrot, y_unrot)

        return lon, lat

    # ------------------------------------------------------------------
    # Grid construction
    # ------------------------------------------------------------------

    @staticmethod
    def _make_edges(centres: np.ndarray, spacing: float) -> np.ndarray:
        """Construct cell edges from uniformly spaced cell centres."""
        edges = np.empty(centres.size + 1, dtype=float)

        edges[1:-1] = 0.5 * (centres[:-1] + centres[1:])

        edges[0] = centres[0] - 0.5 * spacing
        edges[-1] = centres[-1] + 0.5 * spacing

        return edges

    # ------------------------------------------------------------------

    def make_grid(
        self,
        xmin: float,
        xmax: float,
        ymin: float,
        ymax: float,
        dx: float,
        dy: float | None = None,
    ) -> PolarStereoGrid:
        """
        Construct a regular rectangular model grid.

        Parameters
        ----------
        xmin, xmax : float
            Minimum and maximum x coordinate of the OUTER CELL EDGES [m].
        ymin, ymax : float
            Minimum and maximum y coordinate of the OUTER CELL EDGES [m].
        dx, dy : float
            Grid spacing [m].

        Returns
        -------
        PolarStereoGrid
        """
        if dy is None:
            dy = dx

        if dx <= 0:
            msg = "dx must be positive"
            raise ValueError(msg)

        if dy <= 0:
            msg = "dy must be positive"
            raise ValueError(msg)

        if xmax <= xmin:
            msg = "xmax must be greater than xmin"
            raise ValueError(msg)

        if ymax <= ymin:
            msg = "ymax must be greater than ymin"
            raise ValueError(msg)

        # Check that the domain is an integer number of cells.
        nx_float = (xmax - xmin) / dx
        ny_float = (ymax - ymin) / dy

        nx = round(nx_float)
        ny = round(ny_float)

        if not np.isclose(nx_float, nx):
            msg = "(xmax - xmin) must be an integer multiple of dx"
            raise ValueError(msg)

        if not np.isclose(ny_float, ny):
            msg = "(ymax - ymin) must be an integer multiple of dy"
            raise ValueError(msg)

        # Cell centres.
        x = xmin + (np.arange(nx) + 0.5) * dx
        y = ymin + (np.arange(ny) + 0.5) * dy

        # Cell edges.
        xbnds = self._make_edges(x, dx)
        ybnds = self._make_edges(y, dy)

        # 2D cell-centre coordinates.
        X, Y = np.meshgrid(x, y)

        # Geographic coordinates.
        lon, lat = self.xy_to_lonlat(X, Y)

        # Geographic coordinates of the grid corners.
        Xc, Yc = np.meshgrid(xbnds, ybnds)
        lon_corner, lat_corner = self.xy_to_lonlat(Xc, Yc)

        return PolarStereoGrid(
            x=x,
            y=y,
            xbnds=xbnds,
            ybnds=ybnds,
            X=X,
            Y=Y,
            lon=lon,
            lat=lat,
            lon_corner=lon_corner,
            lat_corner=lat_corner,
            dx=dx,
            dy=dy,
        )

    # ------------------------------------------------------------------

    @property
    def proj4(self) -> str:
        """PROJ4 representation of the base CRS."""
        return self.crs.to_proj4()

    # ------------------------------------------------------------------

    @property
    def wkt(self) -> str:
        """WKT representation of the base CRS."""
        return self.crs.to_wkt()
