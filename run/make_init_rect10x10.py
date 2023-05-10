import netCDF4
root = netCDF4.Dataset("init_rect10x10.nc", "w", format="NETCDF4")
metagrp = root.createGroup("structure")
metagrp.type = "devgrid"
datagrp = root.createGroup("data")
xDim = datagrp.createDimension("x", 10)
yDim = datagrp.createDimension("y", 10)
nLay = datagrp.createDimension("nLayers", 1)
hfield_dims = ("y", "x")
mask = datagrp.createVariable("mask", "f8", hfield_dims)
mask[:,::-1] = [[0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,1,1,0,0,0,0],
             [0,0,0,1,1,1,1,0,0,0],
             [0,0,1,1,1,1,1,1,0,0],
             [0,0,1,1,1,1,1,1,0,0],
             [0,0,0,1,1,1,1,1,0,0],
             [0,1,0,1,1,1,1,1,1,0],
             [0,1,1,0,1,1,1,1,1,0],
             [1,1,0,0,1,1,1,0,0,0],
             [1,1,1,1,1,1,0,0,0,0]]
cice = datagrp.createVariable("cice", "f8", hfield_dims)
cice[:,:] = 0.5
hice = datagrp.createVariable("hice", "f8", hfield_dims)
hice[:,:] = 0.1
hsnow = datagrp.createVariable("hsnow", "f8", hfield_dims)
hsnow[:,:] = 0.0
tice = datagrp.createVariable("tice", "f8", ("nLayers", "y", "x",))
tice[:,:,:] = -1.
root.close()
