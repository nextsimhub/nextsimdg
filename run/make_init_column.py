import netCDF4

root = netCDF4.Dataset("init_column.nc", "w", format="NETCDF4")
metagrp = root.createGroup("structure")
metagrp.type = "simple_rectangular"
datagrp = root.createGroup("data")

xDim = datagrp.createDimension("x", 1)
yDim = datagrp.createDimension("y", 1)
nLay = datagrp.createDimension("nLayers", 2)

ice_salinity = 5  # sould match Ice::s in constants.hpp
ocean_salinity = 32.
mu = -0.055  # should match Water::mu in constants.hpp

mask = datagrp.createVariable("mask", "f8", ("x", "y"))
mask[:, :] = [[1]]
cice = datagrp.createVariable("cice", "f8", ("x", "y",))
cice[:, :] = 0.8
hice = datagrp.createVariable("hice", "f8", ("x", "y",))
hice[:, :] = 2
hsnow = datagrp.createVariable("hsnow", "f8", ("x", "y",))
hsnow[:, :] = 0.3
sss = datagrp.createVariable("sss", "f8", ("x", "y",))
sss[:, :] = ocean_salinity
sst = datagrp.createVariable("sst", "f8", ("x", "y",))
sst[:, :] = ocean_salinity * mu
tice = datagrp.createVariable("tice", "f8", ("x", "y", "nLayers"))
tice[:, :, :] = ice_salinity * mu
root.close()
