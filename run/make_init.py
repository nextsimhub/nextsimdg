import netCDF4

# Currently creates the 30 x 30 rect grid files with land-sea mask. For the
# earlier, more rudimentary version of this script, please see the git file
# history. 

nx = 30
ny = 30
nLayers = 1
n_coords = 2

root = netCDF4.Dataset(f"init_rect{nx}x{ny}.nc", "w", format="NETCDF4")

strgrp = root.createGroup("structure")
strgrp.type = "simple_rectangular"

metagrp = root.createGroup("metadata")
metagrp.type = strgrp.type
confgrp = metagrp.createGroup("configuration") # But add nothing to it
timegrp = metagrp.createGroup("time")
time = timegrp.createVariable("time", "i8")
time[:] = 946684800
time.units = "seconds since 1970-01-01T00:00:00Z"
formatted = timegrp.createVariable("formatted", str)
formatted.format = "%Y-%m-%dT%H:%M:%SZ"
formatted[0] = "2000-01-01T00:00:00Z"

datagrp = root.createGroup("data")

x_dim = datagrp.createDimension("xdim", nx)
y_dim = datagrp.createDimension("ydim", ny)
z_dim = datagrp.createDimension("zdim", nLayers)
xvertex_dim = datagrp.createDimension("xvertex", nx + 1)
yvertex_dim = datagrp.createDimension("yvertex", ny + 1)
coords_dim = datagrp.createDimension("ncoords", n_coords)

hfield_dims = ("ydim", "xdim")

mask = datagrp.createVariable("mask", "f8", hfield_dims)
mask[:,::-1] = [[0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0],
             [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0],
             [0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
             [0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
             [0,1,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
             [0,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
             [0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
             [0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
             [0,0,0,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0],
             [0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
             [1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
             [0,1,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
             [0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
             [0,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0],
             [1,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0],
             [1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0],
             [1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
             [1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
             [1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
             [1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0],
             [1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0],
             [1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
             [1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0]]
antimask = 1 - mask[:,:]
cice = datagrp.createVariable("cice", "f8", hfield_dims)
cice[:,::-1] = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,1,1,2,2,2,2,2,2,1,1,2,2,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,1,1,2,3,4,5,4,3,3,2,2,3,2,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,1,2,4,5,6,5,5,4,3,3,4,2,1,1,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,1,2,5,7,8,7,8,6,4,5,3,2,2,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,3,5,7,9,9,9,8,6,8,6,4,2,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,6,0,7,8,9,9,9,9,9,9,8,5,3,2,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,6,0,8,8,9,9,9,9,9,9,9,9,7,5,4,2,1,1,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,7,8,9,9,9,9,9,9,9,9,9,6,5,3,2,1,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,8,9,9,9,9,9,9,9,9,9,9,9,6,5,3,1,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,9,0,8,9,9,9,9,9,9,9,9,9,9,8,0,4,3,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,9,9,9,8,9,9,9,9,9,9,7,5,3,1,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,8,9,9,9,7,6,4,3,2,1,1,1,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,9,9,9,9,9,8,7,5,3,2,2,1,1,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,9,7,5,3,3,4,2,1,1,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,6,5,5,4,3,2,2,2,1,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,4,2,1,1,2,1,1,1,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,3,2,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,3,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
cice[:,:] /= 10
hice = datagrp.createVariable("hice", "f8", hfield_dims)
hice[:,:] = cice[:,:] * 2
hsnow = datagrp.createVariable("hsnow", "f8", hfield_dims)
hsnow[:,:] = cice[:,:] / 2
tice = datagrp.createVariable("tice", "f8", ("zdim", "ydim", "xdim"))
tice[0,:,:] = -0.5 - cice[:,:]

mdi = -3.40282347e38 # Minus float max
# mask data
cice[:,:] = cice[:,:] * mask[:,:] + antimask * mdi
cice.missing_value = mdi
hice[:,:] = hice[:,:] * mask[:,:] + antimask * mdi
hice.missing_value = mdi
hsnow[:,:] = hsnow[:,:] * mask[:,:] + antimask * mdi
hsnow.missing_value = mdi
tice[0,:,:] = tice[0,:,:] * mask[:,:] + antimask * mdi
tice.missing_value = mdi

# coordinates
# element centres
x_var = datagrp.createVariable("x", "f8", hfield_dims)
y_var = datagrp.createVariable("y", "f8", hfield_dims)

d_distance = 150000 # 150 km element spacing

for j in range(0, ny):
    y = (j + 0.5) * d_distance
    for i in range(0, nx):
        x = (i + 0.5) * d_distance
        x_var[j, i] = x
        y_var[j, i] = y

coords = datagrp.createVariable("coords", "f8", ("yvertex", "xvertex", "ncoords"))

for j in range(0, ny + 1):
    y = j * d_distance
    for i in range(0, nx + 1):
        x = i * d_distance
        coords[j, i, 0] = x
        coords[j, i, 1] = y

root.close()
