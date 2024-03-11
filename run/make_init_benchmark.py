import netCDF4
from math import sin
import numpy as np

# Creates initial conditions for the Mehlmann et al. (2021) benchmark case, at 2, 4, 8, and 16 km resolutions.

# Domain size
L = 512
for res in [2, 4, 8, 16]:

    nx = ny = int(L / res)

    # Grid dimensions. Since x and y are switched between the source grid file
    # and the target restart file, the grid dimensions are nfirst and nsecond.
    # nsecond is the size of the dimension that varies fastest.
    nfirst = int(L/res)
    nsecond = int(L/res)
    nLayers = 1
    ncg = 1
    n_dg = 1
    n_dgstress = 3
    n_coords = 2

    fname = f"init_benchmark_{nx}x{ny}.nc"
    print("Producing file", fname)

    root = netCDF4.Dataset(fname, "w", format="NETCDF4")

    structure_name = "parametric_rectangular"
    structgrp = root.createGroup("structure")
    structgrp.type = structure_name

    metagrp = root.createGroup("metadata")
    metagrp.type = structure_name
    confgrp = metagrp.createGroup("configuration") # But add nothing to it
    datagrp = root.createGroup("data")

    nLay = datagrp.createDimension("z", nLayers)
    yDim = datagrp.createDimension("y", nfirst)
    xDim = datagrp.createDimension("x", nsecond)
    yVertexDim = datagrp.createDimension("yvertex", nfirst + 1)
    xVertexDim = datagrp.createDimension("xvertex", nsecond+ 1)
    ycg_dim = datagrp.createDimension("y_cg", nfirst * ncg + 1)
    xcg_dim = datagrp.createDimension("x_cg", nsecond * ncg + 1)
    dg_comp = datagrp.createDimension("dg_comp", n_dg)
    dgs_comp = datagrp.createDimension("dgstress_comp", n_dgstress)
    n_coords_comp = datagrp.createDimension("ncoords", n_coords)

    field_dims = ("y", "x")
    coord_dims = ("yvertex", "xvertex", "ncoords")
    zfield_dims = ("z", "y", "x")

    # Array coordinates
    x = np.zeros((nfirst + 1, nsecond + 1))
    y = np.zeros((nfirst + 1, nsecond + 1))
    '''
    for ix in range(nx):
        x[ix] = (ix+1)*res

    for iy in range(ny):
        y[iy] = (iy+1)*res
    '''
    for ix in range(nfirst + 1):
        for iy in range(nsecond + 1):
            x[ix, iy] = ix*res
            y[ix, iy] = iy*res

    coords = datagrp.createVariable("coords", "f8", coord_dims)
    coords[:,:,0] = x
    coords[:,:,1] = y

    px = np.zeros((nfirst, nsecond))
    py = np.zeros((nfirst, nsecond))
    for ix in range(nfirst):
        for iy in range(nsecond):
            px[ix, iy] = (ix+0.5)*res
            py[ix, iy] = (iy+0.5)*res

    elem_x = datagrp.createVariable("x", "f8", field_dims)
    elem_x[:, :] = px
    elem_y = datagrp.createVariable("y", "f8", field_dims)
    elem_y[:, :] = py

    grid_azimuth = datagrp.createVariable("grid_azimuth", "f8", field_dims)
    grid_azimuth[:, :] = 0.

    # Ice everywhere and all boundaries closed
    mask = datagrp.createVariable("mask", "f8", field_dims)
    mask[:, :] = 1.
    mask[0, :] = 0.
    mask[nx-1, :] = 0.
    mask[:, 0] = 0.
    mask[:, ny-1] = 0.
    antimask = 1 - mask[:, :]

    # Uniform concentration of 100%
    cice = datagrp.createVariable("cice", "f8", field_dims)
    cice[:, :] = 1.

    # Loop over ice thickness to construct the initial conditions. This should be a pattern of slightly varying "hills".
    hice = datagrp.createVariable("hice", "f8", field_dims)

    for ix in range(nx):
        for iy in range(ny):
            hice[ix, iy] = 0.3 + 0.005 * (sin(60e-3 * x[ix, iy]) + sin(30e-3 * y[ix, iy]))

    # Neither snow nor temperature is prescribed. This is a dynamics only benchmark, so temperature is not required.
    hsnow = datagrp.createVariable("hsnow", "f8", field_dims)
    hsnow[:, :] = 0.
    tice = datagrp.createVariable("tice", "f8", ("z", "y", "x"))
    tice[0, :, :] = 0.

    # Ice starts at rest
    u = datagrp.createVariable("u", "f8", field_dims)
    u[:, :] = 0

    v = datagrp.createVariable("v", "f8", field_dims)
    v[:, :] = 0

    # mask data
    mdi = -1. ** 35
    cice[:, :] = cice[:, :] * mask[:, :] + antimask * mdi
    cice.missing_value = mdi
    hice[:, :] = hice[:, :] * mask[:, :] + antimask * mdi
    hice.missing_value = mdi
    hsnow[:, :] = hsnow[:, :] * mask[:, :] + antimask * mdi
    hsnow.missing_value = mdi
    u[:, :] = u[:, :] * mask[:, :] + antimask * mdi
    u.missing_value = mdi
    v[:, :] = v[:, :] * mask[:, :] + antimask * mdi
    v.missing_value = mdi
    tice[0, :, :] = tice[0, :, :] * mask[:, :] + antimask * mdi
    tice.missing_value = mdi

    root.close()
