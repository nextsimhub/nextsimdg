import netCDF4
import numpy as np

# Creates a 128 x 128 ParaGrid restart file with land-sea mask.

nx = 128
ny = 128
nLayers = 3
ncg = 1
n_dg = 1
n_dgstress = 3
n_coords = 2

root = netCDF4.Dataset(f"init_para{nx}x{ny}.nc", "w", format="NETCDF4")

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

hfield_dims = ("y", "x")

mask33 = np.array(
    [[0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
     [0,0,1,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
     [0,0,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0],
     [0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
     [1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
     [1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
     [1,0,1,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
     [0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
     [0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0],
     [1,1,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0],
     [1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
     [1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
     [1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0],
     [1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
     [1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0],
     [1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0],
     [1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
     [1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
     [1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0]])
# scipy isn't working? fine, I'll fake it
mask129x33 = np.zeros((nx+1, 33))
for i in range(33):
    mask129x33[:, i] = np.interp(np.arange(129) / 4, np.arange(33), mask33[:, i])

mask129 = np.zeros((nx+1, ny+1))
for i in range(nx+1):
    mask129[i, :] = np.interp(np.arange(129) / 4, np.arange(33), mask129x33[i, :])

mask = datagrp.createVariable("mask", "f8", hfield_dims)
mask[:,:] = np.rint(mask129[:-1, -2::-1])

antimask = 1 - mask[:,:]

# Some coordinates that don't match the land mask exactly (azimuthal equidistant)
array_size1d = 20.
spacing1d = 2 * array_size1d / nx
limit1d = array_size1d # even number of points + 1
coord1d = np.linspace(-limit1d, limit1d, num=129)

x_coords = np.zeros((nx + 1, ny + 1))
y_coords = np.zeros((nx + 1, ny + 1))
for i in range(nx + 1):
    x_coords[i, :] = coord1d
    y_coords[:, i] = coord1d
    
lat = 90 - (x_coords**2 + y_coords**2)**0.5

lon_x0 = 270. # degrees
lon = np.rad2deg(np.arctan2(y_coords, x_coords)) + lon_x0
# Correct the range to -180 to 180
lon += 180.
lon %= 360.
lon -= 180.

coords = datagrp.createVariable("coords", "f8", ("yvertex", "xvertex", "ncoords"))
coords[:,:,0] = lon
coords[:,:,1] = lat

cice33 = np.array(
    [[0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,1,1,2,2,2,2,2,2,1,1,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,1,2,3,4,5,4,3,3,2,2,3,2,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,2,4,5,6,5,5,4,3,3,4,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,2,5,7,8,7,8,6,4,5,3,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,3,5,7,9,9,9,8,6,8,6,4,2,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,6,0,7,8,9,9,9,9,9,9,8,5,3,2,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,6,0,8,8,9,9,9,9,9,9,9,9,7,5,4,2,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,7,8,9,9,9,9,9,9,9,9,9,6,5,3,2,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,8,9,9,9,9,9,9,9,9,9,9,9,6,5,3,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,9,0,8,9,9,9,9,9,9,9,9,9,9,8,0,4,3,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,9,9,9,8,9,9,9,9,9,9,7,5,3,1,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,8,9,9,9,7,6,4,3,2,1,1,1,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,9,9,9,9,9,8,7,5,3,2,2,1,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,7,5,3,3,4,2,1,1,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,5,5,4,3,2,2,2,1,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,2,1,1,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,2,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,3,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
     )
cice33 /= 10
cice129x33 = np.zeros((nx+1, 33))
for i in range(33):
    cice129x33[:, i] = np.interp(np.arange(129) / 4, np.arange(33), cice33[:, i])

cice129 = np.zeros((nx+1, ny+1))
for i in range(nx+1):
    cice129[i, :] = np.interp(np.arange(129) / 4, np.arange(33), cice129x33[i, :])

cice = datagrp.createVariable("cice", "f8", hfield_dims)
cice[:,:] = cice129[:-1, -2::-1]

hice = datagrp.createVariable("hice", "f8", hfield_dims)
hice[:,:] = cice[:,:] * 2
hsnow = datagrp.createVariable("hsnow", "f8", hfield_dims)
hsnow[:,:] = cice[:,:] / 2
tice = datagrp.createVariable("tice", "f8", ("nLayers", "y", "x"))
tice[0,:,:] = -0.5 - cice[:,:]
tice[1,:,:] = -1.5 - cice[:,:]
tice[2,:,:] = -2.5 - cice[:,:]

mdi = -2.**300
# mask data
cice[:,:] = cice[:,:] * mask[:,:] + antimask * mdi
cice.missing_value = mdi
hice[:,:] = hice[:,:] * mask[:,:] + antimask * mdi
hice.missing_value = mdi
hsnow[:,:] = hsnow[:,:] * mask[:,:] + antimask * mdi
hsnow.missing_value = mdi
for k in range(3):
    tice[k,:,:] = tice[k,:,:] * mask[:,:] + antimask * mdi
tice.missing_value = mdi

root.close()
