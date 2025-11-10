import time
from pathlib import Path

import netCDF4
import numpy as np
from interpolators import topaz4_interpolate
from make_init_base import initMaker


# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(time_struct):
    return f"{topaz_path}/{time_struct.tm_year}/topaz_rean_{time_struct.tm_year}{time_struct.tm_mon:02}.nc"


# Creates a 128 x 128 ParaGrid restart file filled with data from TOPAZ on 2010-01-01
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate an initial state file from TOPAZ4 data and a grid file"
    )
    parser.add_argument(
        "--start-date",
        dest="start_date",
        default="2010-01-01",
        help="Date of the start of the simulation (YYYY-MM-DD)",
    )
    parser.add_argument(
        "--grid-type",
        dest="grid_type",
        default="regular",
        help="The type of input grid to use ('regular' or 'NEMO').",
    )
    parser.add_argument(
        "--grid-file",
        dest="grid_file",
        default="25km_NH.nc",
        help="Path of the input grid file.",
    )
    parser.add_argument(
        "--mask-file",
        dest="mask_file",
        default="25km_NH.nc",
        help="Path of the input mask file.",
    )
    parser.add_argument(
        "--topaz-path",
        dest="topaz_path",
        default=".",
        help="Path containing the TOPAZ4 files.",
    )
    parser.add_argument(
        "--boundary",
        dest="boundary",
        default="open",
        help='One of "open" (default), "closed", or "double".',
    )
    parser.add_argument(
        "--out-suffix",
        dest="out_suffix",
        default="",
        help='Added to the name of the output file before the ending"',
    )

    args = parser.parse_args()
    grid_file = args.grid_file
    topaz_path = args.topaz_path
    boundary = args.boundary
    out_suffix = args.out_suffix

    grid_name = Path(grid_file).stem
    out_name = f"init_{grid_name}{out_suffix}.nc"

    init_base = initMaker(out_name)

    if args.grid_type == "regular":
        init_base.make_geographic_grid(
            grid_name + ".nc", "p", plon_name="plon", plat_name="plat"
        )
        mask = netCDF4.Dataset(grid_name + ".nc", "r")
        init_base.mask[:, :] = mask["mask"][:, :]
    elif args.grid_type == "NEMO":
        init_base.make_geographic_grid(
            grid_name + ".nc",
            "ur",
            plon_name="glamt",
            plat_name="gphit",
            qlon_name="glamf",
            qlat_name="gphif",
        )

        bathy_meter = netCDF4.Dataset(f"{args.mask_file}", "r")
        init_base.mask[:, :] = bathy_meter["Bathymetry"][:, :] > 0.0
    else:
        raise ValueError(f"Grid type {args.grid_type} not supported.")

    if args.boundary in ["closed"]:
        init_base.mask[:, 0] = 0.0
        init_base.mask[:, -1] = 0.0
        init_base.mask[0, :] = 0.0
        init_base.mask[-1, :] = 0.0
    elif args.boundary in ["double"]:
        init_base.mask[:, 0] = 0.0
        init_base.mask[:, -1] = 0.0
        init_base.mask[0, :] = 0.0
        init_base.mask[-1, :] = 0.0
        init_base.mask[:, 1] = 0.0
        init_base.mask[:, -2] = 0.0
        init_base.mask[1, :] = 0.0
        init_base.mask[-2, :] = 0.0

    land_ratio = np.count_nonzero(init_base.mask[:, :]) / init_base.mask[:, :].size
    print(f"ratio of sea (active) cells to total: {land_ratio}")

    # Access the TOPAZ data for initial conditions
    data_time = time.strptime(args.start_date, "%Y-%m-%d")
    source_file_name = topaz4_source_file_name(data_time)
    source_file = netCDF4.Dataset(source_file_name, "r")
    proj_string = source_file["stereographic"].proj4
    source_x = source_file["x"][:]
    source_y = source_file["y"][:]

    element_lon = init_base.get_element_longitude()
    element_lat = init_base.get_element_latitude()

    # Ice concentration and thickness
    init_base.cice[:, :] = topaz4_interpolate(
        element_lon,
        element_lat,
        source_file["siconc"][data_time.tm_mday, :, :].squeeze(),
        source_x,
        source_y,
        proj_string,
    )
    init_base.hice[:, :] = topaz4_interpolate(
        element_lon,
        element_lat,
        source_file["sithick"][data_time.tm_mday, :, :].squeeze(),
        source_x,
        source_y,
        proj_string,
    )

    # Snow thickness
    init_base.hsnow[:, :] = topaz4_interpolate(
        element_lon,
        element_lat,
        source_file["sisnthick"][data_time.tm_mday, :, :].squeeze(),
        source_x,
        source_y,
        proj_string,
    )
    # SSS
    init_base.sss[:, :] = topaz4_interpolate(
        element_lon,
        element_lat,
        source_file["so"][data_time.tm_mday, :, :].squeeze(),
        source_x,
        source_y,
        proj_string,
    )

    # SST
    init_base.sst[:, :] = topaz4_interpolate(
        element_lon,
        element_lat,
        source_file["thetao"][data_time.tm_mday, :, :].squeeze(),
        source_x,
        source_y,
        proj_string,
    )
