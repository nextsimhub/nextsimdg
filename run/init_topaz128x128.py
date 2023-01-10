import netCDF4
import numpy as np
import numpy.ma as ma
import time
import math

topaz_mdi = -32767

# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(field, unix_time):
    unix_tm = time.gmtime(unix_time)
    return f"TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_3m.nc"

# Returns bilinearly interpolated data given array of fractional indices, when some of the data missing
def bilinear_missing(eyes, jays, data, missing):
    i = np.floor(eyes).astype(int)
    j = np.floor(jays).astype(int)
    
    fi = eyes - i
    fj = jays - j

    dataplier = data != missing # False is zero

    weighted_sum = ((1 - fj) * (1 - fi) * data[j, i] * dataplier[j, i] +
        (1 - fj) * (fi) * data[j, i + 1] * dataplier[j, i] +
        (fj) * (1 - fi) * data[j + 1, i] * dataplier[j, i] +
        (fj) * (fi) * data[j + 1, i + 1] * dataplier[j, i])

    sum_of_weights = ((1 - fj) * (1 - fi) * dataplier[j, i] +
        (1 - fj) * (fi) * dataplier[j, i] +
        (fj) * (1 - fi) * dataplier[j, i] +
        (fj) * (fi) * dataplier[j, i])
    
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
    
    nx = 128
    ny = 128
    nLayers = 3
    ncg = 1
    n_dg = 1
    n_dgstress = 3
    n_coords = 2
    
    
    root = netCDF4.Dataset(f"init_topaz{nx}x{ny}.nc", "w", format="NETCDF4")
    
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
    array_size1d = 20. # About 20˚ square
    spacing1d = 2 * array_size1d / nx
    limit1d = array_size1d # even number of points + 1
    coord1d = np.linspace(-limit1d, limit1d, num=129)
    
    x_coords = np.zeros((nx + 1, ny + 1))
    y_coords = np.zeros((nx + 1, ny + 1))
    for i in range(nx + 1):
        x_coords[:, i] = coord1d
        y_coords[i, :] = coord1d
        
    # Polar azimuthal equidistant projection
    # node coordinates
    lat = 90 - (x_coords**2 + y_coords**2)**0.5
    lon = np.rad2deg(np.arctan2(y_coords, x_coords))
    
    # element coordinates
    element_shape = (nx, ny)
    element_lon = np.zeros(element_shape)
    element_lat = np.zeros(element_shape)
    # interpolate lon and lat from nodes to elements, to leave nx x ny arrays
    node_x = np.cos(np.radians(lon)) * np.cos(np.radians(lat))
    node_y = np.sin(np.radians(lon)) * np.cos(np.radians(lat))
    node_z = np.sin(np.radians(lat))
    
    element_x = 0.25 * (node_x[0:-1, 0:-1] + node_x[1:, 0:-1] + node_x[0:-1, 1:] + node_x[1:, 1:])
    element_y = 0.25 * (node_y[0:-1, 0:-1] + node_y[1:, 0:-1] + node_y[0:-1, 1:] + node_y[1:, 1:])
    element_z = 0.25 * (node_z[0:-1, 0:-1] + node_z[1:, 0:-1] + node_z[0:-1, 1:] + node_z[1:, 1:])
    
    element_lon = np.degrees(np.arctan2(element_y, element_x))
    element_lat = np.degrees(np.arctan2(element_z, np.hypot(element_x, element_y)))
    
    # Access the TOPAZ data, initally to get latitudes
    source_file_name = topaz4_source_file_name("hice", data_time)
    source_file = netCDF4.Dataset(topaz4_source_file_name("hice", data_time), "r")
    source_lats = source_file["latitude"][:, :]
    lat_array = source_lats[550:, 380]
    
    # Coordinate values in the file
    nc_lons = datagrp.createVariable("longitude", "f8", ("x", "y"))
    nc_lons[:, :] = element_lon
    nc_lats = datagrp.createVariable("latitude", "f8", ("x", "y"))
    nc_lats[:, :] = element_lat

    coords = datagrp.createVariable("coords", "f8", ("xvertex", "yvertex", "ncoords"))
    coords[:,:,0] = lon
    coords[:,:,1] = lat

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
    
    # Ice temperature
    tice = datagrp.createVariable("tice", "f8", ("x", "y", "z"))
    ice_melt = -0.055 * 5 # Melting point of sea ice (salinity = 5) in ˚C
    # Tice outside the ice pack is the melting point of pure water ice, which is conveniently 0˚C
    ice_temp2d = np.fmin(sst_data, ice_melt) * isice
    tice[:, :, 0] = ice_temp2d
    tice[:, :, 1] = ice_temp2d
    tice[:, :, 2] = ice_temp2d
    