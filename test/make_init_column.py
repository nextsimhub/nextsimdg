import time

import netCDF4
import numpy as np

root = netCDF4.Dataset("init_column.nc", "w", format="NETCDF4")
metagrp = root.createGroup("structure")
structure_name = "parametric_rectangular"
metagrp.type = structure_name
datagrp = root.createGroup("data")

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

nfirst = 1
nsecond = 1
nLayers = 3
ncg = 1
n_dg = 1
n_dgstress = 1
n_coords = 2

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
node_lon = np.zeros((nfirst + 1, nsecond + 1))
node_lat = np.zeros((nfirst + 1, nsecond + 1))

node_lon[0, 0] = 0
node_lon[0, 1] = 270
node_lon[1, 1] = 180
node_lon[1, 0] = 90

lat0 = 89.84  # Gives a box around the pole 25 km a side

node_lat[:, :] = lat0

coords = datagrp.createVariable("coords", "f8", coord_dims)
coords[:, :, 0] = node_lon
coords[:, :, 1] = node_lat

elem_lon = datagrp.createVariable("longitude", "f8", field_dims)
elem_lon[:, :] = 0
elem_lat = datagrp.createVariable("latitude", "f8", field_dims)
elem_lat[:, :] = 90

grid_azimuth = datagrp.createVariable("grid_azimuth", "f8", field_dims)
grid_azimuth[:, :] = 0

ice_salinity = 1  # should match Ice::s in constants.hpp
mu: float = -0.054  # should match Water::mu in constants.hpp
ocean_temperature = -1.89
ocean_salinity = ocean_temperature / mu

mask = datagrp.createVariable("mask", "f8", field_dims)
mask[:, :] = [[1]]
cice = datagrp.createVariable("cice", "f8", field_dims)
cice[:, :] = 1.
hice = datagrp.createVariable("hice", "f8", field_dims)
hice[:, :] = 3.00
hsnow = datagrp.createVariable("hsnow", "f8", field_dims)
hsnow[:, :] = 0.3
sss = datagrp.createVariable("sss", "f8", field_dims)
sss[:, :] = ocean_salinity
sst = datagrp.createVariable("sst", "f8", field_dims)
sst[:, :] = ocean_temperature
tice = datagrp.createVariable("tice", "f8", zfield_dims)
tice[:, :, :] = ice_salinity * mu
# Ice is at rest
u = datagrp.createVariable("u", "f8", field_dims)
u[:, :] = 0
v = datagrp.createVariable("v", "f8", field_dims)
v[:, :] = 0
root.close()
