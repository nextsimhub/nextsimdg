import time

import netCDF4
import numpy as np

# Create a restart file that specifically uses the old dimension names
if __name__ == "__main__":

    # Grid dimensions. x varies fastest
    nfirst = 1
    nsecond = 2
    nLayers = 1
    ncg = 1
    n_dg = 1
    n_dgstress = 1
    n_coords = 2

    root = netCDF4.Dataset("old_names.nc", "w", format="NETCDF4")

    structure_name = "parametric_rectangular"
    structgrp = root.createGroup("structure")
    structgrp.type = structure_name

    metagrp = root.createGroup("metadata")
    metagrp.type = structure_name
    confgrp = metagrp.createGroup("configuration")  # But add nothing to it
    timegrp = metagrp.createGroup("time")
    time_var = timegrp.createVariable("time", "i8")
    data_time = 1263204000
    time_var[:] = data_time
    time.units = "seconds since 1970-01-01T00:00:00Z"
    formatted = timegrp.createVariable("formatted", str)
    formatted.format = "%Y-%m-%dT%H:%M:%SZ"
    formatted[0] = "2010-01-01T00:00:00Z"
    datagrp = root.createGroup("data")

    nLay = datagrp.createDimension("z", nLayers)
    yDim = datagrp.createDimension("y", nfirst)
    xDim = datagrp.createDimension("x", nsecond)
    yVertexDim = datagrp.createDimension("yvertex", nfirst + 1)
    xVertexDim = datagrp.createDimension("xvertex", nsecond + 1)
    ycg_dim = datagrp.createDimension("y_cg", nfirst * ncg + 1)
    xcg_dim = datagrp.createDimension("x_cg", nsecond * ncg + 1)
    dg_comp = datagrp.createDimension("dg_comp", n_dg)
    dgs_comp = datagrp.createDimension("dgstress_comp", n_dgstress)
    n_coords_comp = datagrp.createDimension("ncoords", n_coords)

    field_dims = ("y", "x")
    coord_dims = ("yvertex", "xvertex", "ncoords")
    zfield_dims = ("z", "y", "x")

    # Array coordinates
    space = 25000.0  # 25 km in metres

    x_vertex1d = np.arange(0, space * (nsecond + 1), space)
    y_vertex1d = np.arange(0, space * (nfirst + 1), space)

    x_vertex = np.resize(x_vertex1d, (nfirst + 1, nsecond + 1))
    y_vertex = np.resize(y_vertex1d, (nsecond + 1, nfirst + 1)).transpose()

    coords = datagrp.createVariable("coords", "f8", coord_dims)
    coords[:, :, 0] = x_vertex
    coords[:, :, 1] = y_vertex

    x_cell1d = np.arange(space * 0.5, space * (nsecond + 0.5), space)
    y_cell1d = np.arange(space * 0.5, space * (nfirst + 0.5), space)

    x_cell = np.resize(x_cell1d, (nfirst, nsecond))
    y_cell = np.resize(y_cell1d, (nsecond, nfirst)).transpose()

    elem_x = datagrp.createVariable("x", "f8", field_dims)
    elem_x[:, :] = x_cell
    elem_y = datagrp.createVariable("y", "f8", field_dims)
    elem_y[:, :] = y_cell

    grid_azimuth = datagrp.createVariable("grid_azimuth", "f8", field_dims)
    # Return the grid azimuth to the range -180˚ to 180˚
    grid_azimuth[:, :] = 0.0

    # All sea, everywhere
    mask = datagrp.createVariable("mask", "f8", field_dims)
    mask[:, :] = 1

    # Ice concentration and thickness
    cice = datagrp.createVariable("cice", "f8", field_dims)
    hice = datagrp.createVariable("hice", "f8", field_dims)
    cice[:, :] = 0.96875  # 31/32
    hice[:, :] = 1.5

    # Snow thickness
    hsnow = datagrp.createVariable("hsnow", "f8", field_dims)
    hsnow[:, :] = 0.125

    # SSS
    sss = datagrp.createVariable("sss", "f8", field_dims)
    sss[:, :] = 32.0

    mu = -0.055

    # SST
    sst = datagrp.createVariable("sst", "f8", field_dims)
    sst[:, :] = mu * sss[:, :]

    # Ice temperature
    tice = datagrp.createVariable("tice", "f8", zfield_dims)
    ice_melt = mu * 5  # Melting point of sea ice (salinity = 5) in ˚C
    # Tice outside the ice pack is the melting point of pure water ice, which is conveniently 0˚C
    tice[0, :, :] = ice_melt

    # Ice starts at rest
    u = datagrp.createVariable("u", "f8", field_dims)
    u[:, :] = 0

    v = datagrp.createVariable("v", "f8", field_dims)
    v[:, :] = 0

    root.close()
