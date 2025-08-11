"""
Test the advection of snow as recorded in advection_test.diagnostic.nc.
"""
import os.path
import sys

import netCDF4
import numpy as np

file_name = "advection_test.diagnostic.nc"
# Check that the diagnostic file was output
if not os.path.isfile(file_name):
    print(f"Error: diagnostic file f{file_name} not found.")
    sys.exit(1)

ncFile = netCDF4.Dataset(file_name, "r")

# Check that an hsnow variable exists
hsnow_name = "hsnow"
if hsnow_name not in ncFile.variables:
    print(f"Error: no variable {hsnow_name} found in netCDF path {file_name}.")
    sys.exit(2)

hsnow = ncFile[hsnow_name]

# Check the dimensions of hice
n_hsnow_dim = 4
nx = 50
ny = 21
ndg = 6
if len(hsnow.shape) != n_hsnow_dim:
    print(f"Error: variable {hsnow_name} has the incorrect number of dimensions, expected {n_hsnow_dim}, got {len(hsnow.shape)}.")
    sys.exit(3)
if hsnow.shape[1] != ny:
    print(f"Error: incorrect first dimension for variable {hsnow_name}, expected {ny}, got {hsnow.shape[0]}")
    sys.exit(4)
if hsnow.shape[2] != nx:
    print(f"Error: incorrect second dimension for variable {hsnow_name}, expected {nx}, got {hsnow.shape[1]}")
    sys.exit(5)
if hsnow.shape[3] != ndg:
    print(f"Error: incorrect number of DG components for variable {hsnow_name}, expected {ndg}, got {hsnow.shape[2]}")
    sys.exit(6)

# Extract the x-time plane
snow_plane = hsnow[:, ny//2, :, 0]
# Weight the x coordinate by the snow depth.
x_wgt = np.arange(nx)
wghtd = snow_plane * x_wgt
# Avoid the junk data in the first and last columns to get the weighted mean x position of the snow
x_pos = np.sum(wghtd[:,1:-1], axis = 1) / np.sum(snow_plane[:,1:-1], axis = 1)
# Test that the snow has moved to greater x position over time
if x_pos[0] >= x_pos[-1]:
    print(f"Snow was not advected towards increasing x value")
    sys.exit(7)

print("Success: advection test passed.")
