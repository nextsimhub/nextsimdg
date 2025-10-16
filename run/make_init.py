import netCDF4

# Currently creates the 30 x 30 rect grid files with land-sea mask. For the
# earlier, more rudimentary version of this script, please see the git file
# history.

nx = 30
ny = 30
n_coords = 2

ncFile = netCDF4.Dataset(f"init_rect{nx}x{ny}.nc", "w", format="NETCDF4")

ncFile.structure_name = "simple_rectangular"

time = ncFile.createVariable("time_meta", "i8")
time[:] = 946684800
time.units = "seconds since 1970-01-01T00:00:00Z"
formatted = ncFile.createVariable("formatted", str)
formatted.format = "%Y-%m-%dT%H:%M:%SZ"
formatted[0] = "2000-01-01T00:00:00Z"

x_dim = ncFile.createDimension("xdim", nx)
y_dim = ncFile.createDimension("ydim", ny)
xvertex_dim = ncFile.createDimension("xvertex", nx + 1)
yvertex_dim = ncFile.createDimension("yvertex", ny + 1)
coords_dim = ncFile.createDimension("ncoords", n_coords)

hfield_dims = ("ydim", "xdim")

mask = ncFile.createVariable("mask", "f8", hfield_dims)
# fmt: off
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
# fmt: on
cice[:, :] /= 10
hice = ncFile.createVariable("hice", "f8", hfield_dims)
hice[:, :] = cice[:, :] * 2
hsnow = ncFile.createVariable("hsnow", "f8", hfield_dims)
hsnow[:, :] = cice[:, :] / 2
tsurf = ncFile.createVariable("tsurf", "f8", hfield_dims)
tsurf[:, :] = -0.5 - cice[:, :]

mdi = -3.40282347e38  # Minus float max
# mask data
cice[:, :] = cice[:, :] * mask[:, :] + antimask * mdi
cice.missing_value = mdi
hice[:, :] = hice[:, :] * mask[:, :] + antimask * mdi
hice.missing_value = mdi
hsnow[:, :] = hsnow[:, :] * mask[:, :] + antimask * mdi
hsnow.missing_value = mdi
tsurf[:, :] = tsurf[:, :] * mask[:, :] + antimask * mdi
tsurf.missing_value = mdi

# coordinates
# element centres
x_var = ncFile.createVariable("x", "f8", hfield_dims)
y_var = ncFile.createVariable("y", "f8", hfield_dims)

d_distance = 150000  # 150 km element spacing

for j in range(ny):
    y = (j + 0.5) * d_distance
    for i in range(nx):
        x = (i + 0.5) * d_distance
        x_var[j, i] = x
        y_var[j, i] = y

coords = ncFile.createVariable("coords", "f8", ("yvertex", "xvertex", "ncoords"))

for j in range(ny + 1):
    y = j * d_distance
    for i in range(nx + 1):
        x = i * d_distance
        coords[j, i, 0] = x
        coords[j, i, 1] = y

ncFile.close()
