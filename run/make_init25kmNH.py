import numpy as np
import netCDF4

# Creates a restart file from one of the NH_PS grid files. (Currently hardcoded to 25 km).

grid = netCDF4.Dataset("25km_NH.nc", "r")

nx = grid.dimensions["x"].size
ny = grid.dimensions["y"].size
nLayers = 1
ncg = 1
n_dg = 1
n_dgstress = 1
n_coords = 2

root = netCDF4.Dataset("init_25km_NH.nc", "w", format="NETCDF4")

structure_name = "parametric_rectangular"
structgrp = root.createGroup("structure")
structgrp.type = structure_name

metagrp = root.createGroup("metadata")
metagrp.type = structure_name
confgrp = metagrp.createGroup("configuration") # But add nothing to it
timegrp = metagrp.createGroup("time")
time = timegrp.createVariable("time", "i8")
time[:] = 946684800
time.units = "seconds since 1970-01-01T00:00:00Z"
formatted = timegrp.createVariable("formatted", str)
formatted.format = "%Y-%m-%dT%H:%M:%SZ"
formatted[0] = "2000-01-01T00:00:00Z"

datagrp = root.createGroup("data")

xDim = datagrp.createDimension("x", nx)
yDim = datagrp.createDimension("y", ny)
nLay = datagrp.createDimension("nLayers", nLayers)
xVertexDim = datagrp.createDimension("xvertex", nx + 1)
yVertexDim = datagrp.createDimension("yvertex", ny + 1)
xcg_dim = datagrp.createDimension("x_cg", nx * ncg + 1)
ycg_dim = datagrp.createDimension("y_cg", ny * ncg + 1)
dg_comp = datagrp.createDimension("dg_comp", n_dg)
dgs_comp = datagrp.createDimension("dgstress_comp", n_dgstress)
n_coords_comp = datagrp.createDimension("ncoords", n_coords)

grid_mask = grid["mask"]

mask = datagrp.createVariable("mask", "f8", ("x", "y"))
mask[:,:] = grid_mask[:,:]
antimask = 1 - mask[:,:]

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


cice = datagrp.createVariable("cice", "f8", ("x", "y",))
cice[:,:] = mask[:, :] * 0.95
hice = datagrp.createVariable("hice", "f8", ("x", "y",))
hice[:,:] = cice[:,:] * 2
hsnow = datagrp.createVariable("hsnow", "f8", ("x", "y",))
hsnow[:,:] = cice[:,:] / 2
tice = datagrp.createVariable("tice", "f8", ("x", "y", "nLayers"))
tice[:,:,0] = -0.5 - cice[:,:]

mdi = -2.**300
# mask data
cice[:,:] = cice[:,:] * mask[:,:] + antimask * mdi
cice.missing_value = mdi
hice[:,:] = hice[:,:] * mask[:,:] + antimask * mdi
hice.missing_value = mdi
hsnow[:,:] = hsnow[:,:] * mask[:,:] + antimask * mdi
hsnow.missing_value = mdi
tice[:,:,0] = tice[:,:,0] * mask[:,:] + antimask * mdi
tice.missing_value = mdi

root.close()
