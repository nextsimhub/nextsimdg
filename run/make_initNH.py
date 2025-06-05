import netCDF4
import numpy as np
import numpy.ma as ma
import time
import math
from pathlib import Path

topaz_mdi = -32767

# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(field, unix_time):
    unix_tm = time.gmtime(unix_time)
    if field in ["u", "v"]:
        return f"{topaz_path}/TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_30m.nc"
    else:
        return f"{topaz_path}/TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_3m.nc"

# Returns bilinearly interpolated data given array of fractional indices, when some of the data missing
def bilinear_missing(eyes, jays, data, missing):
    i = np.floor(eyes).astype(int)
    j = np.floor(jays).astype(int)
    
    fi = eyes - i
    fj = jays - j

    dataplier = data != missing # False is zero

    weighted_sum = ((1 - fj) * (1 - fi) * data[j, i] * dataplier[j, i] +
        (1 - fj) * (fi) * data[j, i + 1] * dataplier[j, i + 1] +
        (fj) * (1 - fi) * data[j + 1, i] * dataplier[j + 1, i] +
        (fj) * (fi) * data[j + 1, i + 1] * dataplier[j + 1, i + 1])

    sum_of_weights = ((1 - fj) * (1 - fi) * dataplier[j, i] +
        (1 - fj) * (fi) * dataplier[j, i + 1] +
        (fj) * (1 - fi) * dataplier[j + 1, i] +
        (fj) * (fi) * dataplier[j + 1, i + 1])
    
    weighted_sum += missing * (sum_of_weights == 0)
    sum_of_weights += (sum_of_weights == 0)
    
    return weighted_sum / sum_of_weights


# Returns TOPAZ data interpolated from the data grid and coordinates to the target grid and coordinates
def topaz4_interpolate(target_lon_deg, target_lat_deg, data, lat_array):
    # The TOPAZ grid is assumed and hard coded
    ic = 380
    jc = 550
    
    # Scale of the map and zero longitude
 #   two_r = 1 / math.radians(0.08982849)
    lon0 = math.radians(315.)

#    target_lat = np.radians(target_lat_deg)
    target_lon = np.radians(target_lon_deg)
#    k = two_r * np.cos(target_lat) / np.sqrt(1 + np.sin(target_lat))
    # Use linear interpolation to get the target indices on the topaz grid
    # Negate both latitude arrays so that lat_array is increasing
    topaz_i0 = np.interp(-target_lat_deg, -lat_array, np.arange(len(lat_array)))

    x = topaz_i0 * np.sin(target_lon - lon0)
    y = -topaz_i0 * np.cos(target_lon - lon0)
    target_i = x + ic
    target_j = y + jc
    
    return bilinear_missing(target_i, target_j, data, topaz_mdi)

# Creates a 128 x 128 ParaGrid restart file filled with data from TOPAZ on 2010-01-01
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = "Generate an initial state file from TOPAZ4 data")
    parser.add_argument("--grid-file", dest = "grid_file", default="25km_NH.nc", help = "Path of the NH grid file.")
    parser.add_argument("--topaz-path", dest = "topaz_path", default=".", help = "Path containing the TOPAZ4 files.")
    parser.add_argument("--land-mask", dest = "land_mask", default='data', help='One of "data" or "data_closed_boundary"')
    parser.add_argument("--out-suffix", dest = "out_suffix", default='', help='Added to the name of the output file before the ending"')
    
    args = parser.parse_args()
    grid_file = args.grid_file
    topaz_path = args.topaz_path
    land_mask = args.land_mask
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
    root = netCDF4.Dataset(out_name, "w", format="NETCDF4")
    
    structure_name = "parametric_rectangular"
    structgrp = root.createGroup("structure")
    structgrp.type = structure_name
    
    metagrp = root.createGroup("metadata")
    metagrp.type = structure_name
    confgrp = metagrp.createGroup("configuration") # But add nothing to it
    timegrp = metagrp.createGroup("time")
    time_var = timegrp.createVariable("time", "i8")
    data_time = 1263204000
    time_var[:] = data_time
    time.units = "seconds since 1970-01-01T00:00:00Z"
    formatted = timegrp.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = "2010-01-01T00:00:00Z"
    datagrp = root.createGroup("data")

    yDim = datagrp.createDimension("ydim", nfirst)
    xDim = datagrp.createDimension("xdim", nsecond)
    yVertexDim = datagrp.createDimension("yvertex", nfirst + 1)
    xVertexDim = datagrp.createDimension("xvertex", nsecond+ 1)
    ycg_dim = datagrp.createDimension("y_cg", nfirst * ncg + 1)
    xcg_dim = datagrp.createDimension("x_cg", nsecond * ncg + 1)
    dg_comp = datagrp.createDimension("dg_comp", n_dg)
    dgs_comp = datagrp.createDimension("dgstress_comp", n_dgstress)
    n_coords_comp = datagrp.createDimension("ncoords", n_coords)
    
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
    
    coords = datagrp.createVariable("coords", "f8", coord_dims)
    coords[:,:,0] = node_lon
    coords[:,:,1] = node_lat
    
    elem_lon = datagrp.createVariable("longitude", "f8", field_dims)
    elem_lon[:, :] = grid["plon"][:, :]
    elem_lat = datagrp.createVariable("latitude", "f8", field_dims)
    elem_lat[:, :] = grid["plat"][:, :]
    
    grid_azimuth = datagrp.createVariable("grid_azimuth", "f8", field_dims)
    # Return the grid azimuth to the range -180˚ to 180˚
    grid_azimuth_data = grid["plon"][:, :] + np.degrees(grid["ptheta"][:, :])
    grid_azimuth_data += 180
    grid_azimuth_data %= 360.
    grid_azimuth_data -= 180
    grid_azimuth[:, :] = grid_azimuth_data
    
    # Access the TOPAZ data, initially to get latitudes
    source_file_name = topaz4_source_file_name("hice", data_time)
    source_file = netCDF4.Dataset(source_file_name, "r")
    source_lats = source_file["latitude"][:, :]
    lat_array = source_lats[550:, 380]
    
    # Coordinate values in the file
    element_lon = elem_lon[:, :]
    element_lat = elem_lat[:, :]

    # All fields are stored in one file, already opened as source_file
    # Sea-land mask
    mask = datagrp.createVariable("mask", "f8", field_dims)
    sst_data = topaz4_interpolate(element_lon, element_lat, source_file["temperature"][0, :, :].squeeze(), lat_array)
    mask[:, :] = 1 - ma.getmask(sst_data)
    
    land_ratio = np.count_nonzero(mask) / mask.size
    print(f"ratio of sea (active) cells to total: {land_ratio}")
    if land_mask in ['data_closed_boundary']:
        mask[:,0] = 0.0
        mask[:,-1] = 0.0
        mask[0,:] = 0.0
        mask[-1,:] = 0.0
        land_ratio = np.count_nonzero(mask) / mask.size
        print(f"ratio after adjustment: {land_ratio}")
    
    # Ice concentration and thickness
    cice_data = topaz4_interpolate(element_lon, element_lat, source_file["fice"][0, :, :].squeeze(), lat_array)
    hice_data = topaz4_interpolate(element_lon, element_lat, source_file["hice"][0, :, :].squeeze(), lat_array)

    cice_min = 1e-12
    hice_min = 0.01 # m

    noice = np.logical_or(cice_data < cice_min, hice_data < hice_min)
    isice = 1 - noice
    cice_data *= isice
    
    hice_data *= isice
    hice_data *= cice_data # Convert from ice averaged to grid averaged
    
    cice = datagrp.createVariable("cice", "f8", field_dims)
    hice = datagrp.createVariable("hice", "f8", field_dims)
    cice[:, :] = cice_data
    hice[:, :] = hice_data
    
    # Snow thickness
    hsnow = datagrp.createVariable("hsnow", "f8", field_dims)
    hsnow_data = topaz4_interpolate(element_lon, element_lat, source_file["hsnow"][0, :, :].squeeze(), lat_array)
    hsnow_data *= noice
    hsnow_data *= cice_data
    hsnow[:, :] = hsnow_data

    # SSS
    sss = datagrp.createVariable("sss", "f8", field_dims)
    sss_data = topaz4_interpolate(element_lon, element_lat, source_file["salinity"][0, :, :].squeeze(), lat_array)
    sss[:, :] = sss_data

    mu = -0.055

    # SST
    sst = datagrp.createVariable("sst", "f8", field_dims)
    sst_data = topaz4_interpolate(element_lon, element_lat, source_file["temperature"][0, :, :].squeeze(), lat_array)
    sst[:, :] = sst_data * noice + mu * sss_data * isice

    # Ice temperature
    tsurf = datagrp.createVariable("tsurf", "f8", field_dims)
    tintr = datagrp.createVariable("tinterior", "f8", field_dims)
    tbott = datagrp.createVariable("tbottom", "f8", field_dims)
    #ice_melt = mu * 5 # Melting point of sea ice (salinity = 5) in ˚C
    ice_melt = mu
    # Tice outside the ice pack is the melting point of pure water ice, which is conveniently 0˚C
    ice_temp2d = np.fmin(sst_data, ice_melt)
    tsurf[:, :] = ice_temp2d
    tintr[:, :] = ice_temp2d
    tbott[:, :] = ice_temp2d
    
    uv_source_file = netCDF4.Dataset(topaz4_source_file_name("u", data_time), "r")

    # Ice starts at rest
    u = datagrp.createVariable("u", "f8", field_dims)
    u[:, :] = 0

    v = datagrp.createVariable("v", "f8", field_dims)
    v[:, :] = 0
    
    root.close()
    print(f'Created init file "{out_name}"')
