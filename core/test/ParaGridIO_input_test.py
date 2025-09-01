"""Generate the ParaGridIO_input_test.nc file used by core/test/ParaGrid_test.cpp."""

import netCDF4
import numpy as np

sec_per_hr = 3600
hr_per_day = 24

if __name__ == "__main__":
    nx = 9
    ny = 11
    ncg = 1
    n_dg = 1
    n_dgstress = 1
    n_coords = 2
    element_shape = (nx, ny)
    x1d = np.arange(nx)
    y1d = np.arange(ny)

    element_x = np.tile(x1d, (ny, 1))
    element_y = np.transpose(np.tile(y1d, (nx, 1)))

    time_data = 10 * element_x + element_y
    field_name = "index2d"

    out_file = "ParaGridIO_input_test.nc"
    ncFile = netCDF4.Dataset(out_file, "w", format="NETCDF4")

    ncFile.structure_name = "parametric_rectangular"

    # Use the start time as the timestamp for the file
    formatted = ncFile.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = "1970-01-01T00:00:00Z"
    time_meta = ncFile.createVariable("time_meta", "i8")
    the_time = 0
    time_meta[:] = the_time
    time_meta.units = "seconds since 1970-01-01T00:00:00Z"

    xDim = ncFile.createDimension("x", nx)
    yDim = ncFile.createDimension("y", ny)
    xVertexDim = ncFile.createDimension("xvertex", nx + 1)
    yVertexDim = ncFile.createDimension("yvertex", ny + 1)
    xcg_dim = ncFile.createDimension("x_cg", nx * ncg + 1)
    ycg_dim = ncFile.createDimension("y_cg", ny * ncg + 1)
    dg_comp = ncFile.createDimension("dg_comp", n_dg)
    dgs_comp = ncFile.createDimension("dgstress_comp", n_dgstress)
    n_coords_comp = ncFile.createDimension("ncoords", n_coords)
    tDim = ncFile.createDimension("time", None)

    # Position and time variables
    nc_lons = ncFile.createVariable("longitude", "f8", ("y", "x"))
    nc_lons[:, :] = element_x
    nc_lats = ncFile.createVariable("latitude", "f8", ("y", "x"))
    nc_lats[:, :] = element_y

    nc_times = ncFile.createVariable("time", "f8", ("time"))

    # (unix_times_e, era5_times) = create_era5_times(start_time, stop_time)
    # For each field and time, get the corresponding file name for each dataset
    # get the source data
    data = ncFile.createVariable(field_name, "f8", ("time", "y", "x"))
    data[0, :, :] = time_data
    # 'fill' the time axis
    nc_times[0] = the_time
    ncFile.close()
