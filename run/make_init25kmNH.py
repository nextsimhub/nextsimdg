import netCDF4
import numpy as np
import numpy.ma as ma
import time
import math

import sys # FIXME remove me

topaz_mdi = -32767

# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(field, unix_time):
    unix_tm = time.gmtime(unix_time)
    if field in ["u", "v"]:
        return f"TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_30m.nc"
    else:
        return f"TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_3m.nc"

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
    two_r = 1 / math.radians(0.08982849)
    lon0 = math.radians(315.)

    target_lat = np.radians(target_lat_deg)
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

    grid = netCDF4.Dataset("25km_NH.nc", "r")
    
    nx = grid.dimensions["x"].size
    ny = grid.dimensions["y"].size
    nLayers = 3
    ncg = 1
    n_dg = 1
    n_dgstress = 3
    n_coords = 2
    
    
    root = netCDF4.Dataset("init_25km_NH.nc", "w", format="NETCDF4")
    
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
    
    xDim = datagrp.createDimension("x", nx)
    yDim = datagrp.createDimension("y", ny)
    nLay = datagrp.createDimension("z", nLayers)
    xVertexDim = datagrp.createDimension("xvertex", nx + 1)
    yVertexDim = datagrp.createDimension("yvertex", ny + 1)
    xcg_dim = datagrp.createDimension("x_cg", nx * ncg + 1)
    ycg_dim = datagrp.createDimension("y_cg", ny * ncg + 1)
    dg_comp = datagrp.createDimension("dg_comp", n_dg)
    dgs_comp = datagrp.createDimension("dgstress_comp", n_dgstress)
    n_coords_comp = datagrp.createDimension("ncoords", n_coords)
    
    # Array coordinates
    node_lon = np.zeros((nx + 1, ny + 1))
    node_lat = np.zeros((nx + 1, ny + 1))
    
    node_lon[0:-1, 0:-1] = grid["lon_corners"][:, :, 0]
    node_lon[0:-1, -1] = grid["lon_corners"][:, -1, 1]
    node_lon[-1, -1] = grid["lon_corners"][-1, -1, 2]
    node_lon[-1, 0:-1] = grid["lon_corners"][-1, :, 3]
    
    node_lat[0:-1, 0:-1] = grid["lat_corners"][:, :, 0]
    node_lat[0:-1, -1] = grid["lat_corners"][:, -1, 1]
    node_lat[-1, -1] = grid["lat_corners"][-1, -1, 2]
    node_lat[-1, 0:-1] = grid["lat_corners"][-1, :, 3]
    
    coords = datagrp.createVariable("coords", "f8", ("xvertex", "yvertex", "ncoords"))
    coords[:,:,0] = node_lon
    coords[:,:,1] = node_lat
    
    elem_lon = datagrp.createVariable("longitude", "f8", ("x", "y",))
    elem_lon[:, :] = grid["plon"][:, :]
    elem_lat = datagrp.createVariable("latitude", "f8", ("x", "y",))
    elem_lat[:, :] = grid["plat"][:, :]
    
    grid_azimuth = datagrp.createVariable("grid_azimuth", "f8", ("x", "y"))
    grid_azimuth[:, :] = grid["plon"][:, :] + np.degrees(grid["ptheta"][:, :])
    
    # Access the TOPAZ data, initally to get latitudes
    source_file_name = topaz4_source_file_name("hice", data_time)
    source_file = netCDF4.Dataset(topaz4_source_file_name("hice", data_time), "r")
    source_lats = source_file["latitude"][:, :]
    lat_array = source_lats[550:, 380]
    
    # Coordinate values in the file
    element_lon = elem_lon[:, :]
    element_lat = elem_lat[:, :]

    # All fields are stored in one file, already opened as source_file
    # Sea-land mask
    mask = datagrp.createVariable("mask", "f8", ("x", "y"))
    sst_data = topaz4_interpolate(element_lon, element_lat, source_file["temperature"][0, :, :].squeeze(), lat_array)
    mask[:, :] = 1 - ma.getmask(sst_data)

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
    
    cice = datagrp.createVariable("cice", "f8", ("x", "y"))
    hice = datagrp.createVariable("hice", "f8", ("x", "y"))
    cice[:, :] = cice_data
    hice[:, :] = hice_data
    
    # Snow thickness
    hsnow = datagrp.createVariable("hsnow", "f8", ("x", "y"))
    hsnow_data = topaz4_interpolate(element_lon, element_lat, source_file["hsnow"][0, :, :].squeeze(), lat_array)
    hsnow_data *= noice
    hsnow_data *= cice_data
    hsnow[:, :] = hsnow_data
    
    mu = -0.055
    
    # Ice temperature
    tice = datagrp.createVariable("tice", "f8", ("x", "y", "z"))
    ice_melt = mu * 5 # Melting point of sea ice (salinity = 5) in ˚C
    # Tice outside the ice pack is the melting point of pure water ice, which is conveniently 0˚C
    ice_temp2d = np.fmin(sst_data, ice_melt) * isice
    tice[:, :, 0] = ice_temp2d
    tice[:, :, 1] = ice_temp2d
    tice[:, :, 2] = ice_temp2d
    
    # SSS
    sss = datagrp.createVariable("sss", "f8", ("x", "y"))
    sss_data = topaz4_interpolate(element_lon, element_lat, source_file["salinity"][0, :, :].squeeze(), lat_array)
    sss[:, :] = sss_data

    # SST
    sst = datagrp.createVariable("sst", "f8", ("x", "y"))
    sst_data = topaz4_interpolate(element_lon, element_lat, source_file["temperature"][0, :, :].squeeze(), lat_array)
    sst[:, :] = sst_data * noice + mu * sss_data * isice

    uv_source_file = netCDF4.Dataset(topaz4_source_file_name("u", data_time), "r")

    # Ice starts at rest
    u = datagrp.createVariable("u", "f8", ("x", "y"))
    u_data = topaz4_interpolate(element_lon, element_lat, uv_source_file["u"][0, :, :].squeeze(), lat_array)
    u[:, :] = u_data

    v = datagrp.createVariable("v", "f8", ("x", "y"))
    v_data = topaz4_interpolate(element_lon, element_lat, uv_source_file["v"][0, :, :].squeeze(), lat_array)
    v[:, :] = v_data
    
    root.close()
