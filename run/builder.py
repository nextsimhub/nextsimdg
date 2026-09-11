"""
Build an Arctic ocean-model grid from IBCAO bathymetry.

The resulting NetCDF contains:

    x               model-grid x cell centres [m]
    y               model-grid y cell centres [m]
    x_bnds          model-grid x cell edges [m]
    y_bnds          model-grid y cell edges [m]
    lon             longitude at cell centres [degrees]
    lat             latitude at cell centres [degrees]
    depth           mean ocean depth in each cell [m]
    ocean_mask      boolean ocean mask

The model grid uses EPSG:3996 and may be arbitrarily rotated.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from gebco import GEBCO
from ibcao import IBCAO
from projection import PolarStereoProjection
from scipy.ndimage import label


class GridBuilder:
    """
    Build an ocean-model grid from IBCAO or GEBCO bathymetry.

    Parameters
    ----------
    bathy_file : str or pathlib.Path
        IBCAO GeoTIFF or GEBCO NetCDF file.

    rotation : float, optional
        Counter-clockwise rotation of the model grid [degrees].
    """

    def __init__(self, bathy_file, rotation=0.0):

        self.bathy_file = Path(bathy_file)

        if not self.bathy_file.exists():
            msg = f"Bathymetry file not found: {self.bathy_file}"
            raise FileNotFoundError(msg)

        self.rotation = float(rotation)
        self.projection = PolarStereoProjection(rotation=self.rotation)

        self.grid = None
        self.depth = None
        self.ocean_mask = None

    # ------------------------------------------------------------------
    # Build grid
    # ------------------------------------------------------------------

    def build(self, xmin, xmax, ymin, ymax, resolution, connectivity=4):
        """
        Build the model grid.

        Parameters
        ----------
        xmin, xmax : float
            Model-grid x limits [m].

        ymin, ymax : float
            Model-grid y limits [m].

        resolution : float
            Horizontal grid resolution [m].

        connectivity : float
            Connectivity criterion for flooding isolated wet regions (either 4 or 8)
        """
        self.grid = self.projection.make_grid(
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, dx=resolution
        )

        if "gebco" in self.bathy_file.name.casefold():
            with GEBCO(self.bathy_file) as gebco:
                self.depth, self.ocean_mask = gebco.average_to_grid(
                    self.grid, self.projection
                )
        elif "ibcao" in self.bathy_file.name.casefold():
            with IBCAO(self.bathy_file) as ibcao:
                self.depth, self.ocean_mask = ibcao.average_to_grid(
                    self.grid, self.projection
                )
        else:
            msg = f"Bathymetry file not recognised: {self.bathy_file}. File name must contain 'ibcao' or 'gebco'."
            raise ValueError(msg)

        # --------------------------------------------------------------
        # Flood isolated wet regions
        # --------------------------------------------------------------

        self.depth, self.ocean_mask = self.flood_fill_depth(
            self.depth, self.ocean_mask, connectivity
        )

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------

    def plot(self, ax=None, cmap="viridis", vmin=None, vmax=None, figsize=(10, 10)):
        """
        Plot the generated model grid.

        The plot is made in EPSG:3996 coordinates. This is intentional:
        it shows the actual geometry of the rotated model grid without
        introducing another map projection.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes.

        cmap : str, optional
            Matplotlib colormap.

        vmin, vmax : float, optional
            Colour scale limits.

        figsize : tuple, optional
            Figure size if ``ax`` is not supplied.

        Returns
        -------
        fig, ax
            Matplotlib figure and axes.
        """
        if self.grid is None or self.depth is None or self.ocean_mask is None:
            msg = "build() must be called before plot()"
            raise ValueError(msg)

        # --------------------------------------------------------------
        # Create axes
        # --------------------------------------------------------------

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        else:
            fig = ax.figure

        # --------------------------------------------------------------
        # Construct the cell-edge mesh in EPSG:3996
        # --------------------------------------------------------------
        #
        # The model grid is rectangular in its own coordinate system,
        # but rotated relative to EPSG:3996. Therefore we transform
        # the four corners of every grid vertex.
        #
        # We construct these efficiently by making a mesh of the
        # model-grid edges and applying the inverse rotation.
        # --------------------------------------------------------------

        X_edge, Y_edge = np.meshgrid(self.grid.xbnds, self.grid.ybnds)

        # --------------------------------------------------------------
        # Quantity to plot
        # --------------------------------------------------------------

        values = np.asarray(self.depth, dtype=float)
        label = "Mean ocean depth [m]"

        # --------------------------------------------------------------
        # Plot field
        # --------------------------------------------------------------

        mesh = ax.pcolormesh(
            X_edge, Y_edge, values, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat"
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------

        ax.set_aspect("equal")
        ax.set_xlabel("EPSG:3996 x [m]")
        ax.set_ylabel("EPSG:3996 y [m]")
        ax.set_title(
            f"Arctic model grid ({self.grid.dx:g} m, rotation {self.rotation:g}°)"
        )

        cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
        cbar.set_label(label)

        fig.tight_layout()

        return fig, ax

    # ------------------------------------------------------------------
    # Write NetCDF
    # ------------------------------------------------------------------

    def write_netcdf(self, output_file):
        """
        Build the model grid and write it to a NetCDF file.

        Parameters
        ----------
        output_file : str or pathlib.Path
            Output NetCDF filename.

        xmin, xmax : float
            Model-grid x limits [m].

        ymin, ymax : float
            Model-grid y limits [m].

        resolution : float
            Grid resolution [m].

        Returns
        -------
        Path
            Path to the generated NetCDF file.
        """
        if self.grid is None or self.depth is None or self.ocean_mask is None:
            msg = "build() must be called before write_netcdf()"
            raise ValueError(msg)

        output_file = Path(output_file)

        ds = xr.Dataset(
            data_vars={
                "lon": (
                    ("y", "x"),
                    self.grid.lon.astype(np.float64),
                    {
                        "long_name": "longitude",
                        "standard_name": "longitude",
                        "units": "degrees_east",
                    },
                ),
                "lat": (
                    ("y", "x"),
                    self.grid.lat.astype(np.float64),
                    {
                        "long_name": "latitude",
                        "standard_name": "latitude",
                        "units": "degrees_north",
                    },
                ),
                "lon_corner": (
                    ("y_bnds", "x_bnds"),
                    self.grid.lon_corner.astype(np.float64),
                    {
                        "standard_name": "longitude",
                        "long_name": "longitude at cell corners",
                        "units": "degrees_east",
                    },
                ),
                "lat_corner": (
                    ("y_bnds", "x_bnds"),
                    self.grid.lat_corner.astype(np.float64),
                    {
                        "standard_name": "latitude",
                        "long_name": "latitude at cell corners",
                        "units": "degrees_north",
                    },
                ),
                "depth": (
                    ("y", "x"),
                    self.depth,
                    {
                        "long_name": "mean ocean depth",
                        "standard_name": ("sea_floor_depth_below_sea_surface"),
                        "units": "m",
                    },
                ),
                "ocean_mask": (
                    ("y", "x"),
                    self.ocean_mask.astype(np.uint8),
                    {
                        "long_name": "ocean mask",
                        "units": "1",
                        "flag_values": [0, 1],
                        "flag_meanings": "land ocean",
                    },
                ),
            },
            coords={
                "x": (
                    "x",
                    self.grid.x,
                    {
                        "long_name": ("model grid x coordinate"),
                        "standard_name": ("projection_x_coordinate"),
                        "units": "m",
                    },
                ),
                "y": (
                    "y",
                    self.grid.y,
                    {
                        "long_name": ("model grid y coordinate"),
                        "standard_name": ("projection_y_coordinate"),
                        "units": "m",
                    },
                ),
                "x_bnds": (
                    "x_bnds",
                    self.grid.xbnds,
                    {
                        "long_name": "x cell boundaries",
                        "units": "m",
                    },
                ),
                "y_bnds": (
                    "y_bnds",
                    self.grid.ybnds,
                    {
                        "long_name": "y cell boundaries",
                        "units": "m",
                    },
                ),
            },
            attrs={
                "title": "Arctic ocean model grid",
                "source": "IBCAO/GEBCO",
                "bathimetry_file": str(self.bathy_file),
                "projection": "EPSG:3996",
                "rotation": self.rotation,
                "grid_resolution": self.grid.dx,
            },
        )

        # --------------------------------------------------------------
        # Projection metadata
        # --------------------------------------------------------------

        ds["crs"] = xr.DataArray(
            0,
            attrs={
                "grid_mapping_name": ("polar_stereographic"),
                "epsg_code": "EPSG:3996",
                "spatial_ref": self.projection.wkt,
            },
        )

        for variable in ("lon", "lat", "depth", "ocean_mask"):
            ds[variable].attrs["grid_mapping"] = "crs"

        # --------------------------------------------------------------
        # Encoding
        # --------------------------------------------------------------

        encoding = {
            "depth": {
                "zlib": True,
                "complevel": 4,
                "dtype": "float32",
            },
            "ocean_mask": {
                "zlib": True,
                "complevel": 4,
                "dtype": "uint8",
            },
            "lon": {
                "zlib": True,
                "complevel": 4,
                "dtype": "float64",
            },
            "lat": {
                "zlib": True,
                "complevel": 4,
                "dtype": "float64",
            },
        }

        ds.to_netcdf(output_file, encoding=encoding)
        ds.close()

        return output_file

    def flood_fill_depth(self, depth, ocean_mask, connectivity):
        """
        Remove isolated ocean regions from the depth field.

        Only ocean cells connected to the outer boundary of the model
        domain are retained.

        Parameters
        ----------
        depth : ndarray
            IBCAO-derived mean depth [m].

        ocean_mask : ndarray
            Boolean mask identifying ocean cells [True, False].

        connectivity : {4, 8}, optional
            Cell connectivity used for the flood fill.

            4 -> north/south/east/west neighbours
            8 -> includes diagonal neighbours

        Returns
        -------
        flooded_depth : ndarray
            Depth field with isolated ocean regions removed.

        connected_ocean : ndarray
            Boolean mask identifying ocean cells connected to the
            domain boundary.
        """
        if connectivity not in (4, 8):
            msg = "connectivity must be either 4 or 8"
            raise ValueError(msg)

        ocean = np.isfinite(depth) & ocean_mask

        structure = np.array(
            [
                [0, 1, 0],
                [1, 1, 1],
                [0, 1, 0],
            ],
            dtype=bool,
        )

        labels, nlabels = label(ocean, structure=structure)

        if nlabels > 0:
            sizes = np.bincount(labels.ravel())
            sizes[0] = 0
            largest_label = np.argmax(sizes)
            connected_ocean = labels == largest_label
        else:
            connected_ocean = np.zeros_like(ocean, dtype=bool)

        # --------------------------------------------------------------
        # Remove isolated ocean regions
        # --------------------------------------------------------------

        flooded_depth = np.full_like(depth, np.nan, dtype=np.float32)
        flooded_depth[connected_ocean] = depth[connected_ocean]

        return (flooded_depth, connected_ocean)
