import netCDF4
import numpy as np
import time
import calendar
import math

sec_per_hr = 3600

def create_times(start_tm, stop_tm):
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
    # and TOPAZ4 (hours since 1950-01-01T00:00:00Z Sunday)
    topaz4_epoch = Tm(1950, 1, 1, 0, 0, 0, 6, 1, False)
    topaz4_unix = calendar.timegm(topaz4_epoch)
    topaz4_hours = topaz4_unix / sec_per_hr

    return (unix_times, hour_times - era5_hours, hour_times - topaz4_hours)

def era5_time(unix_time):
    era5_epoch = Tm(1900, 1, 1, 0, 0, 0, 0, 1, False)
    era5_unix = calendar.timegm(era5_epoch)
    era5_sec = unix_time - era5_unix
    return era5_sec / sec_per_hr
    
def era5_source_file_name(field, unix_time):
    file_year = time.gmtime(unix_time).tm_year
    return f"ERA5_{field}_y{file_year}.nc"

def topaz4_source_file_name(field, unix_time):
    unix_tm = time.gmtime(unix_time)
    return f"TP4DAILY_{unix_tm.tm_year}{unix_tm.tm_mon:02}_3m.nc"

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
    (unix_times, era5_times, topaz4_times) = create_times(start_time, stop_time)

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

    atmos_fields = ("dew2m", "lw_in", "sw_in", "pair", "tair", "windspeed")
    era5_fields = ("d2m", "msdwlwrf", "msdwswrf", "msl", "msr", "mtpr", "t2m", "u10", "v10")
    era5_translation = {"dew2m" : "d2m", "lw_in" : "msdwlwrf", "sw_in" : "msdwswrf",
                        "pair" : "msl", "tair" : "t2m"} # windspeed is special
    topaz_fields = ("mlp", "salinity", "temperature")
    topaz_translation = {"mld" : "mlp", "sss" : "salinity", "sst" : "temperature"}

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
    
    # For each field and time, get the corresponding file name for each dataset
    for field_name in atmos_fields:
        data = datagrp.createVariable(field_name, "f8", ("time", "x", "y"))
        if (field_name != "windspeed"):
            era5_field = era5_translation[field_name]
            for target_t_index in range(len(unix_times)):
                # get the source data
                source_file = netCDF4.Dataset(era5_source_file_name(era5_field, unix_times[target_t_index]), "r")
                target_time = era5_times[target_t_index]
                source_times = source_file["time"]
                time_index = target_time - source_times[0]
                source_data = source_file[era5_field][time_index, :, :]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                for i in range(nx):
                    for j in range(ny):
                        time_data[i, j] = 0.
                data[target_t_index, :, :] = time_data
        else:
            for target_t_index in range(len(unix_times)):
                # get the source data
                u_file = netCDF4.Dataset(era5_source_file_name("u10", unix_times[target_t_index]), "r")
                v_file = netCDF4.Dataset(era5_source_file_name("v10", unix_times[target_t_index]), "r")
                target_time = era5_times[target_t_index]
                source_times = u_file["time"]
                time_index = target_time - source_times[0]
                u_data = u_file["u10"][time_index, :, :]
                v_data = v_file["v10"][time_index, :, :]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                for i in range(nx):
                    for j in range(ny):
                        time_data[i, j] = math.hypot(0., 0.)
                data[target_t_index, :, :] = time_data

#    print(era5_source_file_name("d2m", unix_times[0]))
#    print(topaz4_source_file_name("mlp", unix_times[0]))

    era_root.close()