import netCDF4
import numpy as np

# Creates a 128 x 128 ParaGrid-based ERA5 forcings file, but with fake data

nx = 128
ny = 128
nLayers = 3
ncg = 1
n_dg = 1
n_dgstress = 3
n_coords = 2

root = netCDF4.Dataset(f"topaz_test{nx}x{ny}.nc", "w", format="NETCDF4")

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
xVertexDim = datagrp.createDimension("xvertex", nx + 1)
yVertexDim = datagrp.createDimension("yvertex", ny + 1)
n_coords_comp = datagrp.createDimension("ncoords", n_coords)
time_dim = datagrp.createDimension("time", None)

hfield_dims = ("y", "x")
timefield_dims = ("time", "y", "x")

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
    x_coords[:, i] = coord1d
    y_coords[i, :] = coord1d
    
lat = 90 - (x_coords**2 + y_coords**2)**0.5
lon = np.rad2deg(np.arctan2(y_coords, x_coords))

coords = datagrp.createVariable("coords", "f8", ("yvertex", "xvertex", "ncoords"))
coords[:,:,0] = lon
coords[:,:,1] = lat

test_data1d = np.linspace(0, nx - 1, num = nx)
thousandths = np.zeros((nx, ny))
millionths = np.zeros((nx, ny))
for i in range(nx):
    thousandths[i, :] = test_data1d * 1e-6
    millionths[:, i] = test_data1d * 1e-3

test_data = thousandths + millionths
mdi = -2.**300

time_var = datagrp.createVariable("time", "f8", "time")
sst = datagrp.createVariable("sst", "f8", timefield_dims)
sss = datagrp.createVariable("sss", "f8", timefield_dims)
mld = datagrp.createVariable("mld", "f8", timefield_dims)
u = datagrp.createVariable("u", "f8", timefield_dims)
v = datagrp.createVariable("v", "f8", timefield_dims)
ssh = datagrp.createVariable("ssh", "f8", timefield_dims)

# 12 monthly values
for t in range(12):
    time_var[t] = 946684800 + 2592000 * t # 30 day months
    
    sst[t, :, :] = (-test_data - t) * mask + mdi * antimask
    sst.missing_value = mdi
    
    sss[t, :, :] = (30 + test_data + t) * mask + mdi * antimask
    sss.missing_value = mdi
    
    mld[t, :, :] = (test_data + 10 + t) * mask + mdi * antimask
    mld.missing_value = mdi
    
    u[t, :, :] = ((test_data + t) * 0.01) * mask + mdi * antimask
    u.missing_value = mdi
    
    v[t, :, :] = (200 + test_data + 100*t) * mask + mdi * antimask
    v.missing_value = mdi

    ssh[t, :, :] = (test_data + 20 + t) * mask + mdi * antimask
    ssh.missing_value = mdi

root.close()
