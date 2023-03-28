import netCDF4
import numpy as np
import time
import calendar
import math

sec_per_hr = 3600
hr_per_day = 24

zero_C_in_kelvin = 273.15

# Returns arrays of times for the ERA5 file, in Unix and ERA5 format
def create_era5_times(start_tm, stop_tm):
    from collections import namedtuple
    # Define the tm named tuple structure locally
    Tm = namedtuple("Tm", "tm_year tm_mon tm_mday tm_hour tm_min tm_sec tm_wday tm_yday tm_isdst")
    # From tm structures to seconds since Unix epoch
    start_unix = calendar.timegm(start_tm)
    stop_unix = calendar.timegm(stop_tm)
    # Convert to integer hours since epoch
    start_hours = start_unix / sec_per_hr
    stop_hours = stop_unix / sec_per_hr
    hour_times = np.arange(start_hours, stop_hours, 1)
    unix_times = hour_times * sec_per_hr # Yes!
    # Add offsets for ERA5 (hours since 1900-01-01T00:00:00Z Monday)
    era5_epoch = Tm(1900, 1, 1, 0, 0, 0, 0, 1, False)
    era5_unix = calendar.timegm(era5_epoch)
    era5_hours = era5_unix / sec_per_hr

    return (unix_times, hour_times - era5_hours)

# Returns arrays of times for the TOPAZ file, in Unix and TOPAZ format
def create_topaz_times(start_tm, stop_tm):
    from collections import namedtuple
    # Define the tm named tuple structure locally
    Tm = namedtuple("Tm", "tm_year tm_mon tm_mday tm_hour tm_min tm_sec tm_wday tm_yday tm_isdst")
    # From tm structures to seconds since Unix epoch
    start_unix = calendar.timegm(start_tm)
    stop_unix = calendar.timegm(stop_tm)
    # Convert to integer hours since epoch
    start_hours = start_unix / sec_per_hr
    stop_hours = stop_unix / sec_per_hr
    # Topaz only needs one sample every 24 hours
    hour_times = np.arange(start_hours, stop_hours, hr_per_day)
    unix_times = hour_times * sec_per_hr # Yes?
    # Add offsets for TOPAZ (hours since 1950-01-01T00:00:00Z Sunday)
    topaz4_epoch = Tm(1950, 1, 1, 0, 0, 0, 6, 1, False)
    topaz4_unix = calendar.timegm(topaz4_epoch)
    topaz4_hours = topaz4_unix / sec_per_hr

    return (unix_times, hour_times - topaz4_hours)

# Returns the file name that holds the ERA5 data for a given field at a given time
def era5_source_file_name(field, unix_time):
    file_year = time.gmtime(unix_time).tm_year
    return f"ERA5_{field}_y{file_year}.nc"

# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(field, unix_time):
    unix_tm = time.gmtime(unix_time)
    if field in ("u", "v"):
        # Ocean currents come from the 30 m files
        return f"TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_30m.nc"
    else:
        return f"TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_3m.nc"

# Returns bilinearly interpolated data given array of fractional indices
def bilinear(eyes, jays, data):
    i = np.floor(eyes).astype(int)
    j = np.floor(jays).astype(int)
    
    fi = eyes - i
    fj = jays - j

    return ((1 - fj) * (1 - fi) * data[j, i] +
        (1 - fj) * (fi) * data[j, i + 1] +
        (fj) * (1 - fi) * data[j + 1, i] +
        (fj) * (fi) * data[j + 1, i + 1])
    
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

# Returns ERA5 data interpolated from the data grid and coordinates to the target grid and coordinates
def era5_interpolate(target_lons, target_lats, data, data_lons, data_lats):
    target_i = (target_lons - data_lons[0]) / (data_lons[1] - data_lons[0])
    # Make sure that the index is in the range of the size of the longitude array
    target_i %= len(data_lons)

    # Latitudes are on a Gaussian grid, so we need to search a bit.
    target_j = (target_lats - data_lats[0]) / (data_lats[1] - data_lats[0])
    
    return bilinear(target_i, target_j, data)

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
    
    return bilinear_missing(target_i, target_j, data, -32767)
    
# The main script. Calculates ERA5 and TOPAZ forcing files, given a grid to
# interpolate on to and start and stop dates
if __name__ == "__main__":
    # Set up the argument parsing
    import argparse
    parser = argparse.ArgumentParser(description = "Create grid matched forcing files for a ERA5 and TOPAZ4")
    parser.add_argument("--file", dest="file", required = True, help = "A restart file containing the target grid information.")
    parser.add_argument("--start", dest = "start", required = True, help = "The ISO start date for the forcing file.")
    parser.add_argument("--stop", dest = "stop", required = True, help = "The ISO end date for the forcing file.")
    args = parser.parse_args()
    # read the date range
    start_time = time.strptime(args.start, "%Y-%m-%d")
    stop_time = time.strptime(args.stop, "%Y-%m-%d")
    #(unix_times, era5_times, topaz4_times) = create_times(start_time, stop_time)

    # read a grid spec (from a restart file)
    root = netCDF4.Dataset(args.file, "r", format = "NETCDF4")
    structgrp = root.groups["structure"]
    target_structure = "parametric_rectangular"
    if structgrp.type != target_structure:
        print(f"Incorrect structure found: {structgrp.type}, wanted {target_structure}.")
        raise SystemExit
    datagrp = root.groups["data"]
    node_coords = datagrp["coords"]
    # assume lon and lat are 0 and 1 coords
    node_lon = node_coords[:, :, 0]
    node_lat = node_coords[:, :, 1]
    nx = node_lon.shape[0] - 1
    ny = node_lon.shape[1] - 1
    element_shape = (nx, ny)
    element_lon = np.zeros(element_shape)
    element_lat = np.zeros(element_shape)
    # interpolate lon and lat from nodes to elements, to leave nx x ny arrays
    node_x = np.cos(np.radians(node_lon)) * np.cos(np.radians(node_lat))
    node_y = np.sin(np.radians(node_lon)) * np.cos(np.radians(node_lat))
    node_z = np.sin(np.radians(node_lat))
    
    element_x = 0.25 * (node_x[0:-1, 0:-1] + node_x[1:, 0:-1] + node_x[0:-1, 1:] + node_x[1:, 1:])
    element_y = 0.25 * (node_y[0:-1, 0:-1] + node_y[1:, 0:-1] + node_y[0:-1, 1:] + node_y[1:, 1:])
    element_z = 0.25 * (node_z[0:-1, 0:-1] + node_z[1:, 0:-1] + node_z[0:-1, 1:] + node_z[1:, 1:])
    
    element_lon = np.degrees(np.arctan2(element_y, element_x))
    element_lat = np.degrees(np.arctan2(element_z, np.hypot(element_x, element_y)))

    wind_speed = "wind_speed"
    atmos_fields = ("dew2m", "lw_in", "sw_in", "pair", "tair", wind_speed)
    era5_fields = ("d2m", "msdwlwrf", "msdwswrf", "msl", "msr", "mtpr", "t2m", "u10", "v10")
    era5_translation = {"dew2m" : "d2m", "lw_in" : "msdwlwrf", "sw_in" : "msdwswrf",
                        "pair" : "msl", "tair" : "t2m"} # windspeed is special
    kelvin_fields = ("d2m", "t2m")

    ###################################################################

    # ERA5 data
    
    era5_out_file = f"ERA5_{args.start}_{args.stop}.nc"
    era_root = netCDF4.Dataset(era5_out_file, "w", format="NETCDF4")
    structgrp = era_root.createGroup("structure")
    structgrp.type = target_structure
    
    metagrp = era_root.createGroup("metadata")
    metagrp.type = target_structure
    confgrp = metagrp.createGroup("configuration") # But add nothing to it
    timegrp = metagrp.createGroup("time")
    # Use the start time as the timestamp for the file
    formatted = timegrp.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = args.start + "T00:00:00Z"
    time_attr = timegrp.createVariable("time", "i8")
    time_attr[:] = calendar.timegm(start_time)
    time_attr.units = "seconds since 1970-01-01T00:00:00Z"
    
    datagrp = era_root.createGroup("data")
    xDim = datagrp.createDimension("x", nx)
    yDim = datagrp.createDimension("y", ny)
    tDim = datagrp.createDimension("time", None)
    
    # Position and time variables
    nc_lons = datagrp.createVariable("longitude", "f8", ("x", "y"))
    nc_lons[:, :] = element_lon
    nc_lats = datagrp.createVariable("latitude", "f8", ("x", "y"))
    nc_lats[:, :] = element_lat
    
    nc_times = datagrp.createVariable("time", "f8", ("time"))
    
    (unix_times_e, era5_times) = create_era5_times(start_time, stop_time)
    # For each field and time, get the corresponding file name for each dataset
    for field_name in atmos_fields:
        data = datagrp.createVariable(field_name, "f8", ("time", "x", "y"))
        if (field_name != wind_speed):
            era5_field = era5_translation[field_name]
            for target_t_index in range(len(unix_times_e)):
                # get the source data
                source_file = netCDF4.Dataset(era5_source_file_name(era5_field, unix_times_e[target_t_index]), "r")
                source_lons = source_file["longitude"]
                source_lats = source_file["latitude"]
                target_time = era5_times[target_t_index]
                source_times = source_file["time"]
                time_index = target_time - source_times[0]
                source_data = source_file[era5_field][time_index, :, :]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                time_data = era5_interpolate(element_lon, element_lat, source_data, source_lons, source_lats)
                if era5_field in kelvin_fields:
                    time_data -= zero_C_in_kelvin
                data[target_t_index, :, :] = time_data
        else:
            for target_t_index in range(len(unix_times_e)):
                # get the source data
                u_file = netCDF4.Dataset(era5_source_file_name("u10", unix_times_e[target_t_index]), "r")
                v_file = netCDF4.Dataset(era5_source_file_name("v10", unix_times_e[target_t_index]), "r")
                source_lons = u_file["longitude"]
                source_lats = u_file["latitude"]
                target_time = era5_times[target_t_index]
                source_times = u_file["time"]
                time_index = target_time - source_times[0]
                u_data = u_file["u10"][time_index, :, :]
                v_data = v_file["v10"][time_index, :, :]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                time_data = np.hypot(era5_interpolate(element_lon, element_lat, u_data, source_lons, source_lats),
                                     era5_interpolate(element_lon, element_lat, v_data, source_lons, source_lats))
                data[target_t_index, :, :] = time_data
                # Also use the windspeed loop to fill the time axis
                nc_times[time_index] = unix_times_e[target_t_index]
    era_root.close()

    ocean_fields = ("mld", "sss", "sst", "u", "v")
    skip_ocean_fields = ()
    topaz_fields = ("mlp", "salinity", "temperature", "u", "v")
    topaz_translation = {"mld" : "mlp", "sss" : "salinity", "sst" : "temperature", "u" : "u", "v" : "v"}


    ###################################################################
    
    # TOPAZ data
    
    topaz_out_file = f"TOPAZ4_{args.start}_{args.stop}.nc"
    topaz_root = netCDF4.Dataset(topaz_out_file, "w", format="NETCDF4")
    structgrp = topaz_root.createGroup("structure")
    structgrp.type = target_structure
    
    metagrp = topaz_root.createGroup("metadata")
    metagrp.type = target_structure
    confgrp = metagrp.createGroup("configuration") # But add nothing to it
    timegrp = metagrp.createGroup("time")
    # Use the start time as the timestamp for the file
    formatted = timegrp.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = args.start + "T00:00:00Z"
    time_attr = timegrp.createVariable("time", "i8")
    time_attr[:] = calendar.timegm(start_time)
    time_attr.units = "seconds since 1970-01-01T00:00:00Z"
    
    datagrp = topaz_root.createGroup("data")
    xDim = datagrp.createDimension("x", nx)
    yDim = datagrp.createDimension("y", ny)
    tDim = datagrp.createDimension("time", None)
    
    (unix_times_t, topaz4_times) = create_topaz_times(start_time, stop_time)

    source_file = netCDF4.Dataset(topaz4_source_file_name("mlp", unix_times_t[0]), "r")
    source_lats = source_file["latitude"][:, :]
    lat_array = source_lats[550:, 380]
    source_file.close()

    # Position and time variables
    nc_lons = datagrp.createVariable("longitude", "f8", ("x", "y"))
    nc_lons[:, :] = element_lon
    nc_lats = datagrp.createVariable("latitude", "f8", ("x", "y"))
    nc_lats[:, :] = element_lat
    
    nc_times = datagrp.createVariable("time", "f8", ("time"))

    # TOPAZ data is daily, not hourly
    topaz_time_ratio = hr_per_day

    # For each field and time, get the corresponding file name for each dataset
    for field_name in ocean_fields:
        data = datagrp.createVariable(field_name, "f8", ("time", "x", "y"))
        if not field_name in skip_ocean_fields:
            topaz_field = topaz_translation[field_name]
            for target_t_index in range(len(unix_times_t)):
                if field_name == ocean_fields[0]:
                    nc_times[target_t_index] = unix_times_t[target_t_index]
                # get the source data
                source_file = netCDF4.Dataset(topaz4_source_file_name(topaz_field, unix_times_t[target_t_index]), "r")
                target_time = topaz4_times[target_t_index]
                source_times = source_file["time"]
                time_index = (target_time - source_times[0]) // hr_per_day
                source_data = source_file[topaz_field][time_index, :, :].squeeze() # Need to squeeze. Why?
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                time_data = topaz4_interpolate(element_lon, element_lat, source_data, lat_array)
                data[target_t_index, :, :] = time_data
        else:
            for target_t_index in range(len(unix_times_t)):
                # get the source data
                target_time = topaz4_times[target_t_index]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                data[target_t_index, :, :] = time_data
        
    topaz_root.close()
