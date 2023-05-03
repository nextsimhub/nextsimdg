import netCDF4
import numpy as np
import time
import calendar
import math
from errno import ENXIO

sec_per_hr = 3600
hr_per_day = 24

if __name__ == "__main__":
    # Set up the argument parsing
#    import argparse
#    parser = argparse.ArgumentParser(description = "Create fake forcing file for the PAraGridIO test")
#    parser.add_argument("--file", dest="file", required = True, help = "A restart file containing the target grid information.")
#    parser.add_argument("--start", dest = "start", required = True, help = "The ISO start date for the forcing file.")
#    parser.add_argument("--stop", dest = "stop", required = True, help = "The ISO end date for the forcing file.")
#    args = parser.parse_args()
    # read the date range
#    start_time = time.strptime(args.start, "%Y-%m-%d")
#    stop_time = time.strptime(args.stop, "%Y-%m-%d")
    #(unix_times, era5_times, topaz4_times) = create_times(start_time, stop_time)

    # read a grid spec (from a restart file)
#    root = netCDF4.Dataset(args.file, "r", format = "NETCDF4")
#    structgrp = root.groups["structure"]
    target_structure = "parametric_rectangular"
#    if structgrp.type != target_structure:
#        print(f"Incorrect structure found: {structgrp.type}, wanted {target_structure}.")
#        raise SystemExit
#    datagrp = root.groups["data"]
#    node_coords = datagrp["coords"]
    # assume lon and lat are 0 and 1 coords
#    node_lon = node_coords[:, :, 0]
#    node_lat = node_coords[:, :, 1]
    nx = 9
    ny = 11
    element_shape = (nx, ny)
    x1d = np.arange(nx)
    y1d = np.arange(ny)
    
    element_x = np.transpose(np.tile(x1d, (ny, 1)))
    element_y = np.tile(y1d, (nx, 1))
    
    time_data = 10 * element_x + element_y
    field_name = "index2d"
    
    out_file = "ParaGridIO_input_test.nc"
    out_root = netCDF4.Dataset(out_file, "w", format="NETCDF4")
    structgrp = out_root.createGroup("structure")
    structgrp.type = target_structure
    
    metagrp = out_root.createGroup("metadata")
    metagrp.type = target_structure
    confgrp = metagrp.createGroup("configuration") # But add nothing to it
    timegrp = metagrp.createGroup("time")
    # Use the start time as the timestamp for the file
    formatted = timegrp.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = "1970-01-01T00:00:00Z"
    time_attr = timegrp.createVariable("time", "i8")
    the_time = 0
    time_attr[:] = the_time
    time_attr.units = "seconds since 1970-01-01T00:00:00Z"
    
    datagrp = out_root.createGroup("data")
    xDim = datagrp.createDimension("x", nx)
    yDim = datagrp.createDimension("y", ny)
    tDim = datagrp.createDimension("time", None)
    
    # Position and time variables
    nc_lons = datagrp.createVariable("longitude", "f8", ("x", "y"))
    nc_lons[:, :] = element_x
    nc_lats = datagrp.createVariable("latitude", "f8", ("x", "y"))
    nc_lats[:, :] = element_y
    
    nc_times = datagrp.createVariable("time", "f8", ("time"))
    
    #(unix_times_e, era5_times) = create_era5_times(start_time, stop_time)
    # For each field and time, get the corresponding file name for each dataset
        # get the source data
    data = datagrp.createVariable(field_name, "f8", ("time", "x", "y"))
    data[0, :, :] = time_data
    # 'fill' the time axis
    nc_times[0] = the_time
    out_root.close()
