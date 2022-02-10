import netCDF4

nx = 20
ny = 30
nLayers = 3

root = netCDF4.Dataset(f"rect{nx}{ny}.res.nc", "w", format="NETCDF4")

metagrp = root.createGroup("structure")
metagrp.type = "simple_rectangular"


datagrp = root.createGroup("data")

xDim = datagrp.createDimension("x", nx)
yDim = datagrp.createDimension("y", ny)
nLay = datagrp.createDimension("nLayers", nLayers)

cice = datagrp.createVariable("cice", "f8", ("x", "y",))
cice[:,:] = 0.5
hice = datagrp.createVariable("hice", "f8", ("x", "y",))
hice[:,:] = 0.1
hsnow = datagrp.createVariable("hsnow", "f8", ("x", "y",))
hsnow[:,:] = 0.0
sss = datagrp.createVariable("sss", "f8", ("x", "y",))
sss[:,:] = 32.
sst = datagrp.createVariable("sst", "f8", ("x", "y",))
sst[:,:] = -1.
tice = datagrp.createVariable("tice", "f8", ("x", "y", "nLayers"))
tice[:,:,:] = -1.

root.close()
