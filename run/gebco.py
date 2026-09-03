"""
Tiled GEBCO bathymetry handling for a rotated Arctic model grid.

GEBCO is a regular global WGS84 longitude/latitude grid, such as
GEBCO_2026 at 15 arc-second resolution.

The model grid is a regular rectangular grid in EPSG:3996 with
optional rotation in the projected Cartesian plane.

The implementation:

    * reads GEBCO lazily from NetCDF
    * processes the target grid in tiles
    * reads only the required GEBCO region
    * handles arbitrary model-grid rotation
    * handles longitude/dateline crossing
    * handles the North Pole efficiently
    * automatically determines geographic margins
    * displays a progress indicator
    * uses Rasterio's area-weighted Resampling.average

North-pole handling
-------------------

A lon/lat bounding box is particularly inefficient around the pole
because longitude becomes singular there.

Therefore, tiles containing the pole are recursively subdivided.
The recursion continues until the pole-containing tile is sufficiently
small. The remaining non-pole tiles are processed normally.

This prevents a large Arctic model grid from ever requiring the
entire high-latitude GEBCO strip to be loaded into memory.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
from projection import PolarStereoGrid, PolarStereoProjection
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.warp import reproject


class GEBCO:
    """Read and area-average a GEBCO NetCDF file."""

    def __init__(
        self,
        filename,
        tile_size=256,
        pole_tile_size=32,
        margin_factor=1.5,
        progress=True,
    ):
        """
        Initialise the class with the GEBCO filename and some parameters.

        Parameters
        ----------
        filename : str or Path
            GEBCO NetCDF file.

        tile_size : int, optional
            Normal model-grid tile size in cells.

        pole_tile_size : int, optional
            Maximum size of a tile containing the North Pole.
            Such tiles are recursively subdivided until neither
            child contains the pole, or this size is reached.

        margin_factor : float, optional
            Multiplier applied to the automatically calculated
            geographic margin.

        progress : bool, optional
            Display a progress indicator.
        """
        self.filename = Path(filename)

        if not self.filename.exists():
            msg = f"GEBCO file not found: {self.filename}"
            raise FileNotFoundError(msg)

        if tile_size <= 0:
            msg = "tile_size must be positive"
            raise ValueError(msg)

        if pole_tile_size <= 0:
            msg = "pole_tile_size must be positive"
            raise ValueError(msg)

        if pole_tile_size > tile_size:
            msg = "pole_tile_size must not exceed tile_size"
            raise ValueError(msg)

        if margin_factor < 1.0:
            msg = "margin_factor must be >= 1"
            raise ValueError(msg)

        self.tile_size = int(tile_size)
        self.pole_tile_size = int(pole_tile_size)
        self.margin_factor = float(margin_factor)
        self.progress = bool(progress)

        self.ds = xr.open_dataset(
            self.filename,
            engine="netcdf4",
            chunks=None,
            decode_coords="all",
        )

        if "lon" not in self.ds:
            msg = "GEBCO file does not contain a 'lon' coordinate"
            raise ValueError(msg)

        if "lat" not in self.ds:
            msg = "GEBCO file does not contain a 'lat' coordinate"
            raise ValueError(msg)

        if "elevation" not in self.ds:
            msg = "GEBCO file does not contain an 'elevation' variable"
            raise ValueError(msg)

        self.lon = self.ds["lon"]
        self.lat = self.ds["lat"]
        self.elevation = self.ds["elevation"]

        if self.elevation.dims != ("lat", "lon"):
            msg = (
                "Expected elevation dimensions ('lat', 'lon'), "
                f"got {self.elevation.dims}"
            )
            raise ValueError(msg)
        self.crs = "EPSG:4326"

        # Coordinates are tiny compared with elevation and can safely
        # be held in memory.
        self._lon = self.lon.values
        self._lat = self.lat.values

        if self._lon.ndim != 1 or self._lat.ndim != 1:
            msg = "GEBCO longitude and latitude coordinates must be 1-D"
            raise ValueError(msg)

        if self._lon.size < 2 or self._lat.size < 2:
            msg = "GEBCO coordinate arrays are too small"
            raise ValueError(msg)

        if not np.all(np.diff(self._lon) > 0):
            msg = "GEBCO longitude coordinate must be strictly increasing"
            raise ValueError(msg)

        if not np.all(np.diff(self._lat) > 0):
            msg = "GEBCO latitude coordinate must be strictly increasing"
            raise ValueError(msg)

        self.lon_min = float(self._lon[0])
        self.lon_max = float(self._lon[-1])
        self.lat_min = float(self._lat[0])
        self.lat_max = float(self._lat[-1])

        self.dlon = float(np.mean(np.diff(self._lon)))

        self.dlat = float(np.mean(np.diff(self._lat)))

        if not np.allclose(
            np.diff(self._lon),
            self.dlon,
        ):
            msg = "GEBCO longitude grid is not regular"
            raise ValueError(msg)

        if not np.allclose(
            np.diff(self._lat),
            self.dlat,
        ):
            msg = "GEBCO latitude grid is not regular"
            raise ValueError(msg)

    # ------------------------------------------------------------------
    # Resource management
    # ------------------------------------------------------------------

    def close(self):
        self.ds.close()

    def __enter__(self):
        return self

    def __exit__(
        self,
        exc_type,
        exc_value,
        traceback,
    ):
        self.close()

    # ------------------------------------------------------------------
    # Target transform
    # ------------------------------------------------------------------

    @staticmethod
    def _target_transform(
        grid,
        projection,
    ):
        """Construct a Rasterio transform for a model-grid tile."""
        x0, y0 = projection.unrotate(grid.xbnds[0], grid.ybnds[0])

        ex_x, ex_y = projection.unrotate(grid.dx, 0.0)

        ey_x, ey_y = projection.unrotate(0.0, grid.dy)

        x_ul = x0 + grid.ny * ey_x
        y_ul = y0 + grid.ny * ey_y

        return Affine(ex_x, -ey_x, x_ul, ex_y, -ey_y, y_ul)

    # ------------------------------------------------------------------
    # Longitude utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _wrap_lon(lon):
        """Wrap longitude to [-180, 180)."""
        return (np.asarray(lon) + 180.0) % 360.0 - 180.0

    def _gebco_lon(self, lon):
        """Convert longitude to the GEBCO convention."""
        lon = np.asarray(lon)

        if self.lon_min >= 0.0:
            return lon % 360.0

        return self._wrap_lon(lon)

    # ------------------------------------------------------------------
    # Tile perimeter
    # ------------------------------------------------------------------

    def _tile_perimeter(self, grid, i0, i1, j0, j1, n=8):
        """Generate points around the perimeter of a model tile."""
        x0 = grid.xbnds[j0]
        x1 = grid.xbnds[j1]

        y0 = grid.ybnds[i0]
        y1 = grid.ybnds[i1]

        t = np.linspace(0.0, 1.0, n + 1)

        xb = x0 + t * (x1 - x0)
        yb = np.full_like(xb, y0)

        xr = np.full_like(t, x1)
        yr = y0 + t * (y1 - y0)

        xt = x1 - t * (x1 - x0)
        yt = np.full_like(xt, y1)

        xl = np.full_like(t, x0)
        yl = y1 - t * (y1 - y0)

        x = np.concatenate(
            [
                xb,
                xr[1:],
                xt[1:],
                xl[1:],
            ]
        )

        y = np.concatenate(
            [
                yb,
                yr[1:],
                yt[1:],
                yl[1:],
            ]
        )

        return x, y

    def _tile_lonlat(self, grid, projection, i0, i1, j0, j1):
        """Return lon/lat coordinates around a tile."""
        x, y = self._tile_perimeter(grid, i0, i1, j0, j1)

        lon, lat = projection.xy_to_lonlat(x, y)

        return (self._gebco_lon(lon), np.asarray(lat))

    # ------------------------------------------------------------------
    # Pole detection
    # ------------------------------------------------------------------

    def _pole_model_coordinates(self, projection):
        """
        Return the North Pole in model coordinates.

        EPSG:3996 places the pole at (0, 0).
        """
        return projection.rotate(0.0, 0.0)

    def _contains_north_pole(self, grid, projection, i0, i1, j0, j1):
        """Return True if the tile contains the North Pole."""
        px, py = self._pole_model_coordinates(projection)

        return (
            grid.xbnds[j0] <= px <= grid.xbnds[j1]
            and grid.ybnds[i0] <= py <= grid.ybnds[i1]
        )

    # ------------------------------------------------------------------
    # Automatic geographic margin
    # ------------------------------------------------------------------

    def _tile_margin(self, grid, projection, i0, i1, j0, j1):
        """
        Estimate a suitable geographic margin for a tile.

        The minimum margin is two GEBCO pixels. The model-cell scale is
        then converted to degrees, with the longitude scale adjusted
        for latitude.
        """
        _, lat = self._tile_lonlat(grid, projection, i0, i1, j0, j1)

        projected_cell = max(abs(grid.dx), abs(grid.dy))

        metres_per_degree = 111_000.0

        margin_lat = projected_cell / metres_per_degree

        # Longitude becomes singular at the pole. Limit the cosine
        # so that the margin remains finite. Pole tiles are handled
        # separately anyway.
        lat_ref = min(89.9, float(np.min(np.abs(lat))))

        cos_lat = max(np.cos(np.deg2rad(lat_ref)), 1e-3)

        margin_lon = projected_cell / metres_per_degree / cos_lat

        margin = self.margin_factor * max(
            margin_lat, margin_lon, 2.0 * self.dlon, 2.0 * self.dlat
        )

        return margin

    # ------------------------------------------------------------------
    # Longitude intervals
    # ------------------------------------------------------------------

    def _longitude_intervals(self, lon):
        """Return the smallest longitude interval(s) containing lon."""
        lon = self._gebco_lon(np.asarray(lon))

        period = 360.0

        if self.lon_min >= 0.0:
            lon = lon % period
        else:
            lon = self._wrap_lon(lon)

        lon = np.sort(lon)

        extended = np.concatenate(
            [
                lon,
                [lon[0] + period],
            ]
        )

        gaps = np.diff(extended)

        k = int(np.argmax(gaps))

        start = extended[k + 1]
        end = extended[k] + period

        width = end - start

        if width >= period - 1e-10:
            return [
                (
                    self.lon_min,
                    self.lon_max,
                )
            ]

        start %= period
        end %= period

        if self.lon_min < 0.0:
            start = float(self._wrap_lon(start))
            end = float(self._wrap_lon(end))

        else:
            start = float(start)
            end = float(end)

        if start <= end:
            return [(start, end)]

        return [(start, self.lon_max), (self.lon_min, end)]

    # ------------------------------------------------------------------
    # Read GEBCO subset
    # ------------------------------------------------------------------

    def _read_subset(self, lon_intervals, lat_min, lat_max):
        """Read a geographic subset from GEBCO."""
        lat_min = max(self.lat_min, lat_min)

        lat_max = min(self.lat_max, lat_max)

        if lat_max < lat_min:
            msg = "Requested latitude range does not overlap GEBCO"
            raise ValueError(msg)

        arrays = []

        for lon0, lon1 in lon_intervals:
            lon0 = max(self.lon_min, lon0)
            lon1 = min(self.lon_max, lon1)

            if lon1 < lon0:
                continue

            subset = self.elevation.sel(
                lon=slice(lon0, lon1),
                lat=slice(lat_min, lat_max),
            )

            if subset.sizes.get("lon", 0) == 0:
                continue

            arrays.append(subset)

        if not arrays:
            msg = "Requested tile does not overlap GEBCO"
            raise ValueError(msg)

        if len(arrays) == 1:
            return arrays[0]

        # --------------------------------------------------------------
        # Dateline crossing.
        # --------------------------------------------------------------

        first = arrays[0]
        second = arrays[1]

        if second.lon.values[0] < first.lon.values[0]:
            second = second.assign_coords(lon=second.lon + 360.0)
        else:
            first = first.assign_coords(lon=first.lon - 360.0)

        return xr.concat([first, second], dim="lon")

    # ------------------------------------------------------------------
    # Get tile subset
    # ------------------------------------------------------------------

    def _get_tile_subset(
        self,
        grid,
        projection,
        i0,
        i1,
        j0,
        j1,
    ):
        """Determine and read the GEBCO subset for one non-pole tile."""
        lon, lat = self._tile_lonlat(grid, projection, i0, i1, j0, j1)

        margin = self._tile_margin(grid, projection, i0, i1, j0, j1)

        lat_min = max(self.lat_min, float(np.min(lat)) - margin)

        lat_max = min(self.lat_max, float(np.max(lat)) + margin)

        intervals = self._longitude_intervals(lon)

        expanded = []

        for lon0, lon1 in intervals:
            lon0 -= margin
            lon1 += margin

            if self.lon_min < 0.0:
                if lon0 < self.lon_min:
                    expanded.append((lon0 + 360.0, self.lon_max))
                    expanded.append((self.lon_min, lon1))

                elif lon1 > self.lon_max:
                    expanded.append((lon0, self.lon_max))
                    expanded.append((self.lon_min, lon1 - 360.0))

                else:
                    expanded.append((lon0, lon1))

            else:
                lon0 %= 360.0
                lon1 %= 360.0

                if lon0 <= lon1:
                    expanded.append((lon0, lon1))
                else:
                    expanded.append((lon0, 360.0))
                    expanded.append((0.0, lon1))

        return self._read_subset(expanded, lat_min, lat_max)

    # ------------------------------------------------------------------
    # Source transform
    # ------------------------------------------------------------------

    @staticmethod
    def _source_transform(subset):
        """Construct the Rasterio transform for a GEBCO subset."""
        lon = subset.lon.values
        lat = subset.lat.values

        dx = float(np.mean(np.diff(lon)))

        dy = float(np.mean(np.diff(lat)))

        west = lon[0] - 0.5 * dx
        north = lat[-1] + 0.5 * dy

        return Affine(dx, 0.0, west, 0.0, -dy, north)

    # ------------------------------------------------------------------
    # Process one ordinary tile
    # ------------------------------------------------------------------

    def _process_tile(
        self,
        grid,
        projection,
        i0,
        i1,
        j0,
        j1,
    ):
        """Process one tile that does not contain the pole."""

        class Tile:
            pass

        tile = Tile()

        tile.ny = i1 - i0
        tile.nx = j1 - j0

        tile.dx = grid.dx
        tile.dy = grid.dy

        tile.xbnds = grid.xbnds[j0 : j1 + 1]

        tile.ybnds = grid.ybnds[i0 : i1 + 1]

        subset = self._get_tile_subset(grid, projection, i0, i1, j0, j1)

        source = subset.load().values.astype(np.float32, copy=False)

        source = np.flipud(source)

        source_transform = self._source_transform(subset)

        target_transform = self._target_transform(tile, projection)

        elevation = np.full(
            (tile.ny, tile.nx),
            np.nan,
            dtype=np.float32,
        )

        reproject(
            source=source,
            destination=elevation,
            src_transform=source_transform,
            src_crs=self.crs,
            src_nodata=np.nan,
            dst_transform=target_transform,
            dst_crs=projection.crs,
            dst_nodata=np.nan,
            resampling=Resampling.average,
        )

        elevation = np.flipud(elevation)

        return elevation

    # ------------------------------------------------------------------
    # Pole subdivision
    # ------------------------------------------------------------------

    def _split_tile(
        self,
        i0,
        i1,
        j0,
        j1,
    ):
        """Split a tile into up to four approximately equal children."""
        im = (i0 + i1) // 2
        jm = (j0 + j1) // 2

        tiles = []

        # Avoid zero-sized children.
        if i0 < im and j0 < jm:
            tiles.append((i0, im, j0, jm))

        if i0 < im and jm < j1:
            tiles.append((i0, im, jm, j1))

        if im < i1 and j0 < jm:
            tiles.append((im, i1, j0, jm))

        if im < i1 and jm < j1:
            tiles.append((im, i1, jm, j1))

        return tiles

    def _process_tile_recursive(self, grid, projection, i0, i1, j0, j1):
        """
        Process a tile, recursively subdividing it if it contains the North Pole.

        Returns
        -------
        list
            Tuples of:

                (i0, i1, j0, j1, elevation)

        """
        contains_pole = self._contains_north_pole(grid, projection, i0, i1, j0, j1)

        size_y = i1 - i0
        size_x = j1 - j0

        # --------------------------------------------------------------
        # Normal tile.
        # --------------------------------------------------------------

        if not contains_pole:
            elevation = self._process_tile(grid, projection, i0, i1, j0, j1)

            return [(i0, i1, j0, j1, elevation)]

        # --------------------------------------------------------------
        # Pole tile is already small enough.
        #
        # At this point we accept reading all longitudes for this
        # small tile. This is the only place where a full longitude
        # range is deliberately used.
        # --------------------------------------------------------------

        if size_y <= self.pole_tile_size and size_x <= self.pole_tile_size:
            return [self._process_pole_tile(grid, projection, i0, i1, j0, j1)]

        # --------------------------------------------------------------
        # Subdivide.
        # --------------------------------------------------------------

        results = []

        for child in self._split_tile(i0, i1, j0, j1):
            results.extend(self._process_tile_recursive(grid, projection, *child))

        return results

    # ------------------------------------------------------------------
    # Process final small pole tile
    # ------------------------------------------------------------------

    def _process_pole_tile(
        self,
        grid,
        projection,
        i0,
        i1,
        j0,
        j1,
    ):
        """
        Process the final small tile containing the North Pole.

        All longitudes are required in principle, but the latitude
        range is restricted to the tile footprint.

        Because pole tiles are limited by pole_tile_size, this remains
        a bounded-memory operation.
        """

        class Tile:
            pass

        tile = Tile()

        tile.ny = i1 - i0
        tile.nx = j1 - j0

        tile.dx = grid.dx
        tile.dy = grid.dy

        tile.xbnds = grid.xbnds[j0 : j1 + 1]

        tile.ybnds = grid.ybnds[i0 : i1 + 1]

        # Determine the tile latitude range from its perimeter.
        _, lat = self._tile_lonlat(grid, projection, i0, i1, j0, j1)

        margin = self._tile_margin(grid, projection, i0, i1, j0, j1)

        lat_min = max(self.lat_min, float(np.min(lat)) - margin)

        # Because the pole is inside the tile, we need to include
        # the polar point itself.
        lat_max = self.lat_max

        # --------------------------------------------------------------
        # Read only the polar latitude band.
        # --------------------------------------------------------------

        subset = self.elevation.sel(
            lon=slice(self.lon_min, self.lon_max),
            lat=slice(lat_min, lat_max),
        )

        source = subset.load().values.astype(np.float32, copy=False)

        source = np.flipud(source)

        source_transform = self._source_transform(subset)

        target_transform = self._target_transform(tile, projection)

        elevation = np.full((tile.ny, tile.nx), np.nan, dtype=np.float32)

        reproject(
            source=source,
            destination=elevation,
            src_transform=source_transform,
            src_crs=self.crs,
            src_nodata=np.nan,
            dst_transform=target_transform,
            dst_crs=projection.crs,
            dst_nodata=np.nan,
            resampling=Resampling.average,
        )

        elevation = np.flipud(elevation)

        return (i0, i1, j0, j1, elevation)

    # ------------------------------------------------------------------
    # Main operation
    # ------------------------------------------------------------------

    def average_to_grid(self, grid: PolarStereoGrid, projection: PolarStereoProjection):
        """
        Area-average GEBCO onto the model grid.

        The grid is processed in tiles. Tiles containing the North
        Pole are recursively subdivided so that no large geographic
        region surrounding the pole needs to be loaded at once.

        Parameters
        ----------
        grid : PolarStereoGrid
            Target model grid.

        projection : PolarStereoProjection
            Model projection.

        Returns
        -------
        depth : ndarray
            Mean GEBCO elevation [m].
            Ocean values are negative; land is NaN.

        ocean_mask : ndarray
            Boolean ocean mask.
        """
        depth = np.full(grid.shape, np.nan, dtype=np.float32)

        # --------------------------------------------------------------
        # Initial tiles.
        # --------------------------------------------------------------

        initial_tiles = []

        for i0 in range(0, grid.ny, self.tile_size):
            i1 = min(i0 + self.tile_size, grid.ny)

            for j0 in range(0, grid.nx, self.tile_size):
                j1 = min(j0 + self.tile_size, grid.nx)

                initial_tiles.append((i0, i1, j0, j1))

        # --------------------------------------------------------------
        # Determine the number of final tiles.
        #
        # This requires the recursive subdivision to be determined
        # before processing so that the progress indicator has a
        # meaningful total.
        # --------------------------------------------------------------

        final_tiles = []

        def collect(i0, i1, j0, j1):
            contains_pole = self._contains_north_pole(grid, projection, i0, i1, j0, j1)

            size_y = i1 - i0
            size_x = j1 - j0

            if not contains_pole or (
                size_y <= self.pole_tile_size and size_x <= self.pole_tile_size
            ):
                final_tiles.append((i0, i1, j0, j1))
                return

            for child in self._split_tile(i0, i1, j0, j1):
                collect(*child)

        for tile in initial_tiles:
            collect(*tile)

        total_tiles = len(final_tiles)

        # --------------------------------------------------------------
        # Process tiles.
        # --------------------------------------------------------------

        for number, (i0, i1, j0, j1) in enumerate(final_tiles, start=1):
            contains_pole = self._contains_north_pole(grid, projection, i0, i1, j0, j1)

            if contains_pole:
                (_i0, _i1, _j0, _j1, elevation) = self._process_pole_tile(
                    grid, projection, i0, i1, j0, j1
                )

            else:
                elevation = self._process_tile(grid, projection, i0, i1, j0, j1)

            depth[i0:i1, j0:j1] = elevation

            if self.progress:
                percent = 100.0 * number / total_tiles

                print(
                    f"\rGEBCO: tile {number}/{total_tiles} ({percent:5.1f}%)",
                    end="",
                    flush=True,
                )

        if self.progress:
            print()

        ocean_mask = depth < 0.0

        return depth, ocean_mask
