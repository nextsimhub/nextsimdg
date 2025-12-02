import netCDF4
import numpy as np

# Creates a restart file from one of the NH_PS grid files. (Currently hardcoded to 25 km).

grid = netCDF4.Dataset("25km_NH.nc", "r")

nx = grid.dimensions["xdim"].size
ny = grid.dimensions["ydim"].size
ncg = 1
n_dg = 1
n_dgstress = 1
n_coords = 2

ncFile = netCDF4.Dataset("init_25km_NH.nc", "w", format="NETCDF4")

structure_name = "parametric_rectangular"
ncFile.structure_name = structure_name

time = ncFile.createVariable("time_meta", "i8")
time[:] = 946684800
time.units = "seconds since 1970-01-01T00:00:00Z"
formatted = ncFile.createVariable("formatted", str)
formatted.format = "%Y-%m-%dT%H:%M:%SZ"
formatted[0] = "2000-01-01T00:00:00Z"

xDim = ncFile.createDimension("xdim", nx)
yDim = ncFile.createDimension("ydim", ny)
xVertexDim = ncFile.createDimension("xvertex", nx + 1)
yVertexDim = ncFile.createDimension("yvertex", ny + 1)
xcg_dim = ncFile.createDimension("x_cg", nx * ncg + 1)
ycg_dim = ncFile.createDimension("y_cg", ny * ncg + 1)
dg_comp = ncFile.createDimension("dg_comp", n_dg)
dgs_comp = ncFile.createDimension("dgstress_comp", n_dgstress)
n_coords_comp = ncFile.createDimension("ncoords", n_coords)

grid_mask = grid["mask"]

mask = ncFile.createVariable("mask", "f8", ("xdim", "ydim"))
mask[:, :] = grid_mask[:, :]
antimask = 1 - mask[:, :]

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

coords = ncFile.createVariable("coords", "f8", ("xvertex", "yvertex", "ncoords"))
coords[:, :, 0] = node_lon
coords[:, :, 1] = node_lat

elem_lon = ncFile.createVariable(
    "longitude",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
elem_lon[:, :] = grid["plon"][:, :]
elem_lat = ncFile.createVariable(
    "latitude",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
elem_lat[:, :] = grid["plat"][:, :]


cice = ncFile.createVariable(
    "cice",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
cice[:, :] = mask[:, :] * 0.95
hice = ncFile.createVariable(
    "hice",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
hice[:, :] = cice[:, :] * 2
hsnow = ncFile.createVariable(
    "hsnow",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
hsnow[:, :] = cice[:, :] / 2
tsurf = ncFile.createVariable("tsurf", "f8", ("xdim", "ydim"))
tsurf[:, :] = -0.5 - cice[:, :]
tbott = ncFile.createVariable("tbottom", "f8", ("xdim", "ydim"))
tbott[:, :] = -1.8
tintr = ncFile.createVariable("tinterior", "f8", ("xdim", "ydim"))
tintr[:, :] = 0.5 * (tsurf[:, :] + tbott[:, :])

sst = ncFile.createVariable(
    "sst",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
sst[:, :] = -cice[:, :]
sss = ncFile.createVariable(
    "sss",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
sss[:, :] = cice[:, :] * 33.68
u = ncFile.createVariable(
    "u",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
u[:, :] = 0.0
v = ncFile.createVariable(
    "v",
    "f8",
    (
        "xdim",
        "ydim",
    ),
)
v[:, :] = 0.0

# velocity in the middle of the domain
midx = (nx + 1) // 2
midy = (nx + 1) // 2

for i in range(nx + 1):
    for j in range(ny + 1):
        if (np.abs(j - midy) < 20) and (np.abs(i - midx) < 20):
            u[i, j] = 0.1
            v[i, j] = 0.1


"""
mdi =  -2.**300
# mask data
cice[:,:] = cice[:,:] * mask[:,:] + antimask * mdi
cice.missing_value = mdi
hice[:,:] = hice[:,:] * mask[:,:] + antimask * mdi
hice.missing_value = mdi
hsnow[:,:] = hsnow[:,:] * mask[:,:] + antimask * mdi
hsnow.missing_value = mdi
tice[:,:,0] = tice[:,:,0] * mask[:,:] + antimask * mdi
tice.missing_value = mdi
sst[:,:] = sst[:,:] * mask[:,:] + antimask * mdi
sst.missing_value = mdi
sss[:,:] = sss[:,:] * mask[:,:] + antimask * mdi
sss.missing_value = mdi
u[:,:] = u[:,:] * mask[:,:] + antimask * mdi
u.missing_value = mdi
v[:,:] = v[:,:] * mask[:,:] + antimask * mdi
v.missing_value = mdi
"""
ncFile.close()
