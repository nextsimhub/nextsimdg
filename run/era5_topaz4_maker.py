import argparse
import calendar
import time
from collections import namedtuple

import netCDF4
import numpy as np
from interpolators import (
    era5_interpolate,
    heading_to_greenland,
    rotate_velocities,
    topaz4_interpolate,
)

sec_per_hr = 3600
hr_per_day = 24

zero_C_in_kelvin = 273.15


def create_era5_times(start_tm, stop_tm):
    """
    Calculate times for the ERA5 file.

    :param start_tm: Start time for the time series
    :param stop_tm: End time for the time series
    :return: Arrays of times for the ERA5 file, in Unix and ERA5 format
    """
    # Define the tm named tuple structure locally
    Tm = namedtuple(
        "Tm", "tm_year tm_mon tm_mday tm_hour tm_min tm_sec tm_wday tm_yday tm_isdst"
    )
    # From tm structures to seconds since Unix epoch
    start_unix = calendar.timegm(start_tm)
    stop_unix = calendar.timegm(stop_tm)
    # Convert to integer hours since epoch
    start_hours = start_unix / sec_per_hr
    stop_hours = stop_unix / sec_per_hr
    hour_times = np.arange(start_hours, stop_hours, 1)
    unix_times = hour_times * sec_per_hr  # Yes!
    # Add offsets for ERA5 (hours since 1900-01-01T00:00:00Z Monday)
    era5_epoch = Tm(1900, 1, 1, 0, 0, 0, 0, 1, False)
    era5_unix = calendar.timegm(era5_epoch)
    era5_hours = era5_unix / sec_per_hr

    return (unix_times, hour_times - era5_hours)


def create_topaz_times(start_tm, stop_tm):
    """
    Calculate times for the TOPAZ file.

    :param start_tm: Start time for the time series
    :param stop_tm: End time for the time series
    :return: Arrays of times for the TOPAZ file, in Unix and TOPAZ format
    """
    # Define the tm named tuple structure locally
    Tm = namedtuple(
        "Tm", "tm_year tm_mon tm_mday tm_hour tm_min tm_sec tm_wday tm_yday tm_isdst"
    )
    # From tm structures to seconds since Unix epoch
    start_unix = calendar.timegm(start_tm)
    stop_unix = calendar.timegm(stop_tm)
    # Convert to integer hours since epoch
    start_hours = start_unix / sec_per_hr
    stop_hours = stop_unix / sec_per_hr
    # Topaz only needs one sample every 24 hours
    hour_times = np.arange(start_hours, stop_hours, hr_per_day)
    unix_times = hour_times * sec_per_hr  # Yes?
    # Add offsets for TOPAZ (hours since 1950-01-01T00:00:00Z Sunday)
    topaz4_epoch = Tm(1950, 1, 1, 0, 0, 0, 6, 1, False)
    topaz4_unix = calendar.timegm(topaz4_epoch)
    topaz4_hours = topaz4_unix / sec_per_hr

    return (unix_times, hour_times - topaz4_hours)


def era5_source_file_name(field, unix_time, path):
    """
    Construct the file name that holds the ERA5 data for a given field at a given time.

    :param field: Name of the variable to be read in from the ERA5 file
    :param unix_time: Time in Unix format
    :param path: Path for the ERA5 file
    :return: File name
    """
    file_year = time.gmtime(unix_time).tm_year
    return f"{path}/ERA5_{field}_y{file_year}.nc"


# Returns the file name that holds the TOPAZ data for a given field at a given time
def topaz4_source_file_name(unix_time, path):
    """
    Construct the file name that holds the TOPAZ4 data for a given field at a given time.

    :param field: Name of the variable to be read in from the TOPAZ4 file
    :param unix_time: Time in Unix format
    :param path: Path for the TOPAZ4 file
    :return: File name
    """
    unix_tm = time.gmtime(unix_time)
    return (
        f"{path}/{unix_tm.tm_year}/topaz_rean_{unix_tm.tm_year}{unix_tm.tm_mon:02}.nc"
    )


if __name__ == "__main__":
    """
    The main script. Calculates ERA5 and TOPAZ forcing files, given a grid to interpolate on to and start and stop dates.
    """

    # Set up the argument parsing
    parser = argparse.ArgumentParser(
        description="Create grid matched forcing files for a ERA5 and TOPAZ4"
    )
    parser.add_argument(
        "--file",
        dest="file",
        required=True,
        help="A restart file containing the target grid information.",
    )
    parser.add_argument(
        "--start",
        dest="start",
        required=True,
        help="The ISO start date for the forcing file.",
    )
    parser.add_argument(
        "--stop",
        dest="stop",
        required=True,
        help="The ISO end date for the forcing file.",
    )
    parser.add_argument(
        "--prefix",
        dest="prefix",
        required=False,
        help="A string to prefix the created files with.",
    )
    parser.add_argument(
        "--era5-path",
        dest="era5_path",
        default="",
        help="Path for the ERA5 forcing files",
    )
    parser.add_argument(
        "--topaz4-path",
        dest="topaz4_path",
        default="",
        help="Path for the TOPAZ4 forcing files",
    )
    args = parser.parse_args()
    # read the date range
    start_time = time.strptime(args.start, "%Y-%m-%d")
    stop_time = time.strptime(args.stop, "%Y-%m-%d")
    # (unix_times, era5_times, topaz4_times) = create_times(start_time, stop_time)

    if args.prefix is not None:
        filepfx = args.prefix + "."
    else:
        filepfx = ""

    era5_path = args.era5_path
    topaz4_path = args.topaz4_path

    # read a grid spec (from a restart file)
    ncFile = netCDF4.Dataset(args.file, "r", format="NETCDF4")

    target_structure = "parametric_rectangular"
    if ncFile.structure_name != target_structure:
        sys_err = (
            f"Incorrect structure found: {ncFile.structure_name}, wanted "
            f"{target_structure}."
        )
        raise SystemExit(sys_err)
    node_coords = ncFile["coords"]
    # assume lon and lat are 0 and 1 coords
    node_lon = node_coords[:, :, 0]
    node_lat = node_coords[:, :, 1]
    ny = node_lon.shape[0] - 1
    nx = node_lon.shape[1] - 1
    element_shape = (ny, nx)
    element_lon = np.zeros(element_shape)
    element_lat = np.zeros(element_shape)
    # interpolate lon and lat from nodes to elements, to leave ny x nx arrays
    node_x = np.cos(np.radians(node_lon)) * np.cos(np.radians(node_lat))
    node_y = np.sin(np.radians(node_lon)) * np.cos(np.radians(node_lat))
    node_z = np.sin(np.radians(node_lat))

    element_x = 0.25 * (
        node_x[0:-1, 0:-1] + node_x[1:, 0:-1] + node_x[0:-1, 1:] + node_x[1:, 1:]
    )
    element_y = 0.25 * (
        node_y[0:-1, 0:-1] + node_y[1:, 0:-1] + node_y[0:-1, 1:] + node_y[1:, 1:]
    )
    element_z = 0.25 * (
        node_z[0:-1, 0:-1] + node_z[1:, 0:-1] + node_z[0:-1, 1:] + node_z[1:, 1:]
    )

    # take grid coordinates out of the init file
    element_lon[:, :] = ncFile["longitude"][:, :]
    element_lat[:, :] = ncFile["latitude"][:, :]
    # manual computation is slightly off and can create inconsistent forcings
    #    element_lon = np.degrees(np.arctan2(element_y, element_x))
    #    element_lat = np.degrees(np.arctan2(element_z, np.hypot(element_x, element_y)))

    # azimuth of the +y grid direction in radians
    element_azimuth = np.radians(element_lon * 0 + 45)

    wind_speed = "wind_speed"
    atmos_fields = ("dew2m", "lw_in", "sw_in", "pair", "tair", wind_speed)
    era5_fields = (
        "d2m",
        "msdwlwrf",
        "msdwswrf",
        "msl",
        "msr",
        "mtpr",
        "t2m",
        "u10",
        "v10",
    )
    era5_translation = {
        "dew2m": "d2m",
        "lw_in": "msdwlwrf",
        "sw_in": "msdwswrf",
        "pair": "msl",
        "tair": "t2m",
    }  # windspeed is special
    kelvin_fields = ("d2m", "t2m")

    ###################################################################

    # ERA5 data

    era5_out_file = f"{filepfx}ERA5_{args.start}_{args.stop}.nc"
    print(f"Writing ERA5 data to {era5_out_file}")
    era5_ncFile = netCDF4.Dataset(era5_out_file, "w", format="NETCDF4")
    era5_ncFile.structure_name = target_structure

    xDim = era5_ncFile.createDimension("x_dim", nx)
    yDim = era5_ncFile.createDimension("y_dim", ny)
    tDim = era5_ncFile.createDimension("time", None)

    hfield_dims = ("y_dim", "x_dim")
    timefield_dims = ("time", "y_dim", "x_dim")

    # Position and time variables
    nc_lons = era5_ncFile.createVariable("longitude", "f8", hfield_dims)
    nc_lons[:, :] = element_lon
    nc_lats = era5_ncFile.createVariable("latitude", "f8", hfield_dims)
    nc_lats[:, :] = element_lat

    greenland_headings = heading_to_greenland(element_lat, element_lon)

    nc_times = era5_ncFile.createVariable("time", "f8", ("time"))
    nc_times.units = "seconds since 1970-01-01T00:00:00Z"

    (unix_times_e, era5_times) = create_era5_times(start_time, stop_time)
    # For each field and time, get the corresponding file name for each dataset
    for field_name in atmos_fields:
        data = era5_ncFile.createVariable(field_name, "f8", timefield_dims)
        if field_name != wind_speed:
            era5_field = era5_translation[field_name]
            for target_t_index in range(len(unix_times_e)):
                # get the source data
                source_file = netCDF4.Dataset(
                    era5_source_file_name(
                        era5_field, unix_times_e[target_t_index], era5_path
                    ),
                    "r",
                )
                source_lons = source_file["longitude"]
                source_lats = source_file["latitude"]
                target_time = era5_times[target_t_index]
                source_times = source_file["time"]
                time_index = target_time - source_times[0]
                source_data = source_file[era5_field][time_index, :, :]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((ny, nx))
                time_data = era5_interpolate(
                    element_lon, element_lat, source_data, source_lons, source_lats
                )
                if era5_field in kelvin_fields:
                    time_data -= zero_C_in_kelvin
                data[target_t_index, :, :] = time_data
        else:
            # Also handle the wind components along with the wind speed
            u_var = era5_ncFile.createVariable("u", "f8", timefield_dims)
            v_var = era5_ncFile.createVariable("v", "f8", timefield_dims)
            for target_t_index in range(len(unix_times_e)):
                # get the source data
                u_file = netCDF4.Dataset(
                    era5_source_file_name(
                        "u10", unix_times_e[target_t_index], era5_path
                    ),
                    "r",
                )
                v_file = netCDF4.Dataset(
                    era5_source_file_name(
                        "v10", unix_times_e[target_t_index], era5_path
                    ),
                    "r",
                )
                source_lons = u_file["longitude"]
                source_lats = u_file["latitude"]
                target_time = era5_times[target_t_index]
                source_times = u_file["time"]
                time_index = target_time - source_times[0]
                u_data_source = u_file["u10"][time_index, :, :]
                v_data_source = v_file["v10"][time_index, :, :]
                # Fix the polar row winds by copying the pole-adjacent row
                # Pole in the first row
                if np.isclose(abs(source_lats[0]), 90.0, atol=0.125):
                    u_data_source[0, :] = u_data_source[1, :]
                    v_data_source[0, :] = u_data_source[1, :]
                # Pole in the last row
                if np.isclose(abs(source_lats[-1]), 90.0, atol=0.125):
                    u_data_source[-1, :] = u_data_source[-2, :]
                    v_data_source[-1, :] = u_data_source[-2, :]
                # Now interpolate the source data to the target grid
                u_data_target = np.zeros((ny, nx))
                u_data_target = era5_interpolate(
                    element_lon, element_lat, u_data_source, source_lons, source_lats
                )
                v_data_target = np.zeros((ny, nx))
                v_data_target = era5_interpolate(
                    element_lon, element_lat, v_data_source, source_lons, source_lats
                )
                speed_data = np.hypot(u_data_target, v_data_target)
                data[target_t_index, :, :] = speed_data
                # Rotate the components from the ERA5 geographic grid to the Greenland displaced pole grid
                (u_data, v_data) = rotate_velocities(
                    u_data_target, -v_data_target, greenland_headings
                )
                u_var[target_t_index, :, :] = u_data
                v_var[target_t_index, :, :] = -v_data
                # Also use the windspeed loop to fill the time axis
                nc_times[time_index] = unix_times_e[target_t_index]
    era5_ncFile.close()

    ocean_fields = ("mld", "sss", "sst", "ssh")
    skip_ocean_fields = ()
    topaz_translation = {"mld": "mlotst", "sss": "so", "sst": "thetao", "ssh": "zos"}

    ###################################################################

    # TOPAZ data

    topaz_out_file = f"{filepfx}TOPAZ4_{args.start}_{args.stop}.nc"
    topaz_ncFile = netCDF4.Dataset(topaz_out_file, "w", format="NETCDF4")
    print(f"Writing TOPAZ4 data to {topaz_out_file}")
    topaz_ncFile.structure_name = target_structure

    xDim = topaz_ncFile.createDimension("x_dim", nx)
    yDim = topaz_ncFile.createDimension("y_dim", ny)
    tDim = topaz_ncFile.createDimension("time", None)

    (unix_times_t, topaz4_times) = create_topaz_times(start_time, stop_time)

    source_file = netCDF4.Dataset(
        topaz4_source_file_name(unix_times_t[0], topaz4_path), "r"
    )
    source_lats = source_file["latitude"][:, :]
    lat_array = source_lats[550:, 380]
    source_file.close()

    # Position and time variables
    nc_lons = topaz_ncFile.createVariable("longitude", "f8", hfield_dims)
    nc_lons[:, :] = element_lon
    nc_lats = topaz_ncFile.createVariable("latitude", "f8", hfield_dims)
    nc_lats[:, :] = element_lat

    nc_times = topaz_ncFile.createVariable("time", "f8", ("time"))
    nc_times.units = "seconds since 1970-01-01T00:00:00Z"

    # TOPAZ data is daily, not hourly
    topaz_time_ratio = hr_per_day

    # The current components are offset by 45˚
    topaz_phi0 = 45  # degrees

    # For each field and time, get the corresponding file name for each dataset
    for field_name in ocean_fields:
        data = topaz_ncFile.createVariable(field_name, "f8", timefield_dims)
        if field_name not in skip_ocean_fields:
            topaz_field = topaz_translation[field_name]
            for target_t_index in range(len(unix_times_t)):
                if field_name == ocean_fields[0]:
                    nc_times[target_t_index] = unix_times_t[target_t_index]
                # get the source data
                source_file = netCDF4.Dataset(
                    topaz4_source_file_name(unix_times_t[target_t_index], topaz4_path),
                    "r",
                )
                target_time = topaz4_times[target_t_index]
                source_times = source_file["time"]
                time_index = (target_time - source_times[0]) // hr_per_day
                # Index the time and squeeze the time dimension away
                source_data = source_file[topaz_field][time_index, :, :].squeeze()
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                proj_string = source_file["stereographic"].proj4
                source_x = source_file["x_dim"][:]
                source_y = source_file["y_dim"][:]
                time_data = topaz4_interpolate(
                    element_lon,
                    element_lat,
                    source_data,
                    source_x,
                    source_y,
                    proj_string,
                )
                data[target_t_index, :, :] = time_data
        else:
            for target_t_index in range(len(unix_times_t)):
                # get the source data
                target_time = topaz4_times[target_t_index]
                # Now interpolate the source data to the target grid
                time_data = np.zeros((nx, ny))
                data[target_t_index, :, :] = time_data

    # Ocean currents
    udata = topaz_ncFile.createVariable("u", "f8", timefield_dims)
    vdata = topaz_ncFile.createVariable("v", "f8", timefield_dims)
    for target_t_index in range(len(unix_times_t)):
        u_source_file = netCDF4.Dataset(
            topaz4_source_file_name(unix_times_t[target_t_index], topaz4_path), "r"
        )
        v_source_file = netCDF4.Dataset(
            topaz4_source_file_name(unix_times_t[target_t_index], topaz4_path), "r"
        )
        target_time = topaz4_times[target_t_index]
        source_times = u_source_file["time"]
        time_index = (target_time - source_times[0]) // hr_per_day
        u_source_data = u_source_file["vxo"][
            time_index, :, :
        ].squeeze()  # Need to squeeze. Why?
        v_source_data = v_source_file["vyo"][time_index, :, :].squeeze()
        # We want the ocean velocity to be zero on land - not some interpolated value
        u_source_data[np.isnan(u_source_data)] = 0.0
        v_source_data[np.isnan(v_source_data)] = 0.0
        u_source_data_tgrid = np.zeros((nx, ny))
        v_source_data_tgrid = np.zeros((nx, ny))
        # Interpolate the current components on the TOPAZ basis on to the new grid
        proj_string = source_file["stereographic"].proj4
        source_x = source_file["x_dim"][:]
        source_y = source_file["y_dim"][:]
        u_source_data_tgrid = topaz4_interpolate(
            element_lon, element_lat, u_source_data, source_x, source_y, proj_string
        )
        v_source_data_tgrid = topaz4_interpolate(
            element_lon, element_lat, v_source_data, source_x, source_y, proj_string
        )

        # Rotate from grid coordinates to geographic eastward/northward components
        rotation_rad = np.radians(element_lon + topaz_phi0 - element_azimuth)
        rotation_rad += heading_to_greenland(element_lat, element_lon)

        (u_target_data, v_target_data) = rotate_velocities(
            u_source_data_tgrid, v_source_data_tgrid, rotation_rad
        )

        udata[target_t_index, :, :] = u_target_data
        vdata[target_t_index, :, :] = v_target_data

    topaz_ncFile.close()
