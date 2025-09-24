import math
import time

import netCDF4
import numpy as np
from numpy import ma

ice_dim = 16
ice_start = 2
snow_offset = 4
snow_hw = 4
ice_stop = ice_start + ice_dim + 1

def get_data(name):
    nxNH = 154
    nyNH = 121
    nx = 50
    ny = 21
    x0 = 30
    y0 = 30

    lon0 = 0.
    dlon = 0.25
    lat0 = -2.5
    dlat = 0.25


    ncoords = 2
    if name == "nx":
        return nx
    elif name == "ny":
        return ny
    elif name == "ncoords":
        return ncoords
    elif name == "longitude":
        lon1d = np.arange(nx) * dlon + lon0
        c_base = np.ones((ny, nx))
        return c_base * lon1d
    elif name == "latitude":
        lat1d = np.arange(ny) * dlat + lat0
        c_base = np.ones((ny, nx))
        return c_base * lat1d[:, np.newaxis]
    elif name == "coords":
        vlon1d = (np.arange(nx + 1) - 0.5) * dlon + lon0
        vlat1d = (np.arange(ny + 1) - 0.5) * dlat + lat0
        v_base = np.ones((ny + 1, nx + 1))
        coords = np.ndarray((ny + 1, nx + 1, 2))
        coords[:, :, 0] = v_base * vlon1d
        coords[:, :, 1] = v_base * vlat1d[:, np.newaxis]
        return coords
    elif name == "grid_azimuth":
        return np.zeros((ny, nx))
    elif name == "mask":
        mask_data = np.ones((ny, nx))
        mask_data[0, :] = 0
        mask_data[-1, :] = 0
        mask_data[:, 0] = 0
        mask_data[:, -1] = 0
        return mask_data
    elif name == "cice":
        cice = np.ones((ny, nx))
        cice[:ice_start, :] = 0.
        cice[:, :ice_start] = 0.
        cice[ice_start + ice_dim + 1:, :] = 0.
        cice[:, ice_start + ice_dim + 1:] = 0.
        return cice
    elif name =="hice":
        return get_data("cice")
    elif name == "sss":
        return np.ones((ny, nx)) * 35
    elif name =="sst":
        return np.ones((ny, nx)) * 35 * -0.0555
    elif name == "tsurf":
        return get_data("cice") * -1.5
if __name__ == "__main__":

    # Grid dimensions
    nx = get_data("nx")
    ny = get_data("ny")
    ncg = 1
    n_dg = 1
    n_dgstress = 3
    n_coords = get_data("ncoords")

    ncFile = netCDF4.Dataset("advection_test_init.nc", "w", format="NETCDF4")

    structure_name = "parametric_rectangular"
    ncFile.structure_name = structure_name

    time_var = ncFile.createVariable("time_meta", "i8")
    data_time = 1263204000
    time_var[:] = data_time
    time.units = "seconds since 1970-01-01T00:00:00Z"
    formatted = ncFile.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = "2010-01-01T00:00:00Z"

    yDim = ncFile.createDimension("ydim", ny)
    xDim = ncFile.createDimension("xdim", nx)
    yVertexDim = ncFile.createDimension("yvertex", ny + 1)
    xVertexDim = ncFile.createDimension("xvertex", nx + 1)
    ycg_dim = ncFile.createDimension("y_cg", ny * ncg + 1)
    xcg_dim = ncFile.createDimension("x_cg", nx * ncg + 1)
    dg_comp = ncFile.createDimension("dg_comp", n_dg)
    dgs_comp = ncFile.createDimension("dgstress_comp", n_dgstress)
    n_coords_comp = ncFile.createDimension("ncoords", n_coords)

    field_dims = ("ydim", "xdim")
    coord_dims = ("yvertex", "xvertex", "ncoords")

    ncFile.createVariable("coords", "f8", coord_dims)[:] = get_data("coords")
    ncFile.createVariable("longitude", "f8", field_dims)[:] = get_data("longitude")
    ncFile.createVariable("latitude", "f8", field_dims)[:] = get_data("latitude")

    ncFile.createVariable("grid_azimuth", "f8", field_dims)[:] = get_data("grid_azimuth")

    ncFile.createVariable("mask", "f8", field_dims)[:, :] = get_data("mask")

    ncFile.createVariable("cice", "f8", field_dims)[:, :] = get_data("cice")
    ncFile.createVariable("hice", "f8", field_dims)[:, :] = get_data("hice")

    # Snow thickness. Zero everywhere on the ocean
    hsnow = 0 * get_data("hice")
    snow_start = ice_start + snow_offset
    snow_stop = ice_stop - snow_offset
    hsnow[snow_start:snow_stop, snow_start:snow_stop] = 0.75
    ncFile.createVariable("hsnow", "f8", field_dims)[:, :] = hsnow

    # SSS
    ncFile.createVariable("sss", "f8", field_dims)[:, :] = get_data("sss")

    # SST
    ncFile.createVariable("sst", "f8", field_dims)[:, :] = get_data("sst")

    # Ice temperature
    ncFile.createVariable("tsurf", "f8", field_dims)[:, :] = get_data("tsurf")

    # Ice starts at rest
    ncFile.createVariable("u", "f8", field_dims)[:, :] = 0

    ncFile.createVariable("v", "f8", field_dims)[:, :] = 0

    ncFile.close()
