import netCDF4
root = netCDF4.Dataset("dev1.res.nc", "w", format="NETCDF4")
metagrp = root.createGroup("structure")
metagrp.type = "devgrid"
datagrp = root.createGroup("data")
xDim = datagrp.createDimension("x", 10)
yDim = datagrp.createDimension("y", 10)
nLay = datagrp.createDimension("nLayers", 1)
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
