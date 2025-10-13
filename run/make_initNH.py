import time
from pathlib import Path

import netCDF4
import numpy as np
from interpolators import topaz4_interpolate
from numpy import ma

topaz_mdi = -32767

# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(unix_time):
    unix_tm = time.gmtime(unix_time)
    return f"{topaz_path}/{unix_tm.tm_year}/topaz_rean_{unix_tm.tm_year}{unix_tm.tm_mon:02}.nc"

# Creates a 128 x 128 ParaGrid restart file filled with data from TOPAZ on 2010-01-01
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = "Generate an initial state file from TOPAZ4 data")
    parser.add_argument("--grid-file", dest = "grid_file", default="25km_NH.nc", help = "Path of the NH grid file.")
    parser.add_argument("--topaz-path", dest = "topaz_path", default=".", help = "Path containing the TOPAZ4 files.")
    parser.add_argument("--boundary", dest = "boundary", default="open", help='One of "open" (default) or "closed"')
    parser.add_argument("--out-suffix", dest = "out_suffix", default="", help='Added to the name of the output file before the ending"')

    args = parser.parse_args()
    grid_file = args.grid_file
    topaz_path = args.topaz_path
    boundary = args.boundary
    out_suffix = args.out_suffix

    grid = netCDF4.Dataset(f"{grid_file}", "r")

    # Grid dimensions. Since x and y are switched between the source grid file
    # and the target restart file, the grid dimensions are nfirst and nsecond.
    # nsecond is the size of the dimension that varies fastest.
    nfirst = grid.dimensions["x"].size
    nsecond = grid.dimensions["y"].size
    print(f"grid size: {nfirst} x {nsecond}")
    ncg = 1
    n_dg = 1
    n_dgstress = 3
    n_coords = 2

    grid_name = Path(grid_file).stem
    out_name = f"init_{grid_name}{out_suffix}.nc"
    ncFile = netCDF4.Dataset(out_name, "w", format="NETCDF4")

    structure_name = "parametric_rectangular"
    ncFile.structure_name = structure_name

    time_var = ncFile.createVariable("time_meta", "i8")
    data_time = 1263204000
    time_var[:] = data_time
    time.units = "seconds since 1970-01-01T00:00:00Z"
    formatted = ncFile.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = "2010-01-01T00:00:00Z"

    yDim = ncFile.createDimension("ydim", nfirst)
    xDim = ncFile.createDimension("xdim", nsecond)
    yVertexDim = ncFile.createDimension("yvertex", nfirst + 1)
    xVertexDim = ncFile.createDimension("xvertex", nsecond+ 1)
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

    node_lon[0:-1, 0:-1] = grid["lon_corners"][:, :, 0]
    node_lon[0:-1, -1] = grid["lon_corners"][:, -1, 1]
    node_lon[-1, -1] = grid["lon_corners"][-1, -1, 2]
    node_lon[-1, 0:-1] = grid["lon_corners"][-1, :, 3]

    node_lat[0:-1, 0:-1] = grid["lat_corners"][:, :, 0]
    node_lat[0:-1, -1] = grid["lat_corners"][:, -1, 1]
    node_lat[-1, -1] = grid["lat_corners"][-1, -1, 2]
    node_lat[-1, 0:-1] = grid["lat_corners"][-1, :, 3]

    coords = ncFile.createVariable("coords", "f8", coord_dims)
    coords[:,:,0] = node_lon
    coords[:,:,1] = node_lat

    elem_lon = ncFile.createVariable("longitude", "f8", field_dims)
    elem_lon[:, :] = grid["plon"][:, :]
    elem_lat = ncFile.createVariable("latitude", "f8", field_dims)
    elem_lat[:, :] = grid["plat"][:, :]

    grid_azimuth = ncFile.createVariable("grid_azimuth", "f8", field_dims)
    # Return the grid azimuth to the range -180˚ to 180˚
    grid_azimuth_data = grid["plon"][:, :] + np.degrees(grid["ptheta"][:, :])
    grid_azimuth_data += 180
    grid_azimuth_data %= 360.
    grid_azimuth_data -= 180
    grid_azimuth[:, :] = grid_azimuth_data

    # Access the TOPAZ data, initially to get coordinates and projection
    source_file_name = topaz4_source_file_name(data_time)
    source_file = netCDF4.Dataset(source_file_name, "r")
    proj_string = source_file["stereographic"].proj4
    source_x = source_file["x"][:]
    source_y = source_file["y"][:]

    # Coordinate values in the file
    element_lon = elem_lon[:, :]
    element_lat = elem_lat[:, :]

    # All fields are stored in one file, already opened as source_file
    # Sea-land mask
    mask = ncFile.createVariable("mask", "f8", field_dims)
    mask[:,:] = grid["mask"][:,:]

    land_ratio = np.count_nonzero(mask) / mask.size
    print(f"ratio of sea (active) cells to total: {land_ratio}")
    if boundary in ["closed"]:
        mask[:,0] = 0.0
        mask[:,-1] = 0.0
        mask[0,:] = 0.0
        mask[-1,:] = 0.0
        land_ratio = np.count_nonzero(mask) / mask.size
        print(f"ratio after adjustment: {land_ratio}")

    nanmask = np.where(mask[:,:] == 0, np.nan, 1)

    # Ice concentration and thickness
    cice_data = topaz4_interpolate(element_lon, element_lat, source_file["siconc"][0, :, :].squeeze(), source_x,
                                   source_y, proj_string)
    hice_data = topaz4_interpolate(element_lon, element_lat, source_file["sithick"][0, :, :].squeeze(), source_x,
                                   source_y, proj_string)

    cice_min = 1e-12
    hice_min = 0.01 # m

    noice = np.logical_or(cice_data < cice_min, hice_data < hice_min)
    isice = 1 - noice
    cice_data *= isice

    hice_data *= isice
    hice_data *= cice_data # Convert from ice averaged to grid averaged

    cice = ncFile.createVariable("cice", "f8", field_dims)
    hice = ncFile.createVariable("hice", "f8", field_dims)
    cice[:, :] = nanmask * cice_data
    hice[:, :] = nanmask * hice_data

    # Snow thickness
    hsnow = ncFile.createVariable("hsnow", "f8", field_dims)
    hsnow_data = topaz4_interpolate(element_lon, element_lat, source_file["sisnthick"][0, :, :].squeeze(), source_x,
                                    source_y, proj_string)
    hsnow_data *= noice
    hsnow_data *= cice_data
    hsnow[:, :] = nanmask * hsnow_data

    # SSS
    sss = ncFile.createVariable("sss", "f8", field_dims)
    sss_data = topaz4_interpolate(element_lon, element_lat, source_file["so"][0, :, :].squeeze(), source_x, source_y,
                                  proj_string)
    sss[:, :] = nanmask * sss_data

    mu = -0.055

    # SST
    sst = ncFile.createVariable("sst", "f8", field_dims)
    sst_data = topaz4_interpolate(element_lon, element_lat, source_file["thetao"][0, :, :].squeeze(), source_x,
                                  source_y, proj_string)
    sst[:, :] = nanmask * sst_data * noice + mu * sss_data * isice

    # Ice starts at rest
    u = ncFile.createVariable("u", "f8", field_dims)
    u[:, :] = 0

    v = ncFile.createVariable("v", "f8", field_dims)
    v[:, :] = 0

    ncFile.close()
    print(f'Created init file "{out_name}"')
