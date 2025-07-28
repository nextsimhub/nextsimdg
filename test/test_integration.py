import numpy as np
import netCDF4
import sys
import os.path

file_name = "out.integration_test.nc"
# Check that the restart file was output
if not os.path.isfile(file_name):
    print(f"Error: restart file f{file_name} not found.")
    sys.exit(1)

root = netCDF4.Dataset(file_name, "r")
data_group_name = "data"

# Check the data group exists
if not data_group_name in root.groups:
    print(f"Error: no group {data_group_name} found in file {file_name}.")
    sys.exit(2)

data_group = root[data_group_name]

# Check that an hice variable exists
hice_name = "hice"
if not hice_name in data_group.variables:
    print(f"Error: no varaible {hice_name} found in netCDF path {file_name}/{data_group_name}.")
    sys.exit(3)

hice = data_group[hice_name]

# Check the dimensions of hice
n_hice_dim = 3
nx = 154
ny = 121
ndg = 6
if len(hice.shape) != n_hice_dim:
    print(f"Error: variable {hice_name} has the incorrect number of dimensions, expected {n_hice_dim}, got {len(hice.shape)}.")
    sys.exit(4)
if hice.shape[0] != ny:
    print(f"Error: incorrect first dimension for variable {hice_name}, expected {ny}, got {hice.shape[0]}")
    sys.exit(5)
if hice.shape[1] != nx:
    print(f"Error: incorrect second dimension for variable {hice_name}, expected {nx}, got {hice.shape[1]}")
    sys.exit(6)
if hice.shape[2] != ndg:
    print(f"Error: incorrect number of DG components for variable {hice_name}, expected {ndg}, got {hice.shape[2]}")
    sys.exit(7)
    
# Test the contents of hice
x_lo = 51
x_hi = x_lo*2 + 1
y_lo = 40
y_hi = y_lo*2 + 1
dg0 = 0
dg1 = 1

# DG0, average or finite difference values
hice_dg0 = hice[y_lo:y_hi:y_lo, x_lo:x_hi:x_lo, dg0]
print(f"hice_dg0={hice_dg0}")
# These values should be greater than zero and less than 10 m
hice_min = 0
hice_max = 10
if np.any(hice_dg0 <= hice_min) or np.any(hice_dg0 > hice_max):
    print(f"Error: ice thickness at at least one sample point is outside the range 0 < h_ice < 10 m.")
    sys.exit(8)
hice_dg1 = hice[y_lo:y_hi:y_lo, x_lo:x_hi:x_lo, dg1]
abs_dg_min = 0
abs_dg_max = 1
# Only check if all the DG1 components are zero. Some are allowed to be 0. But any above the maximum is a failure
if np.all(np.abs(hice_dg1) <= abs_dg_min) or np.any(np.abs(hice_dg1) > abs_dg_max):
    print(f"Error: ice thickness DG component {dg1} out of acceptable bounds.")
    sys.exit(9)
    
cice_name = "cice"
cice = data_group[cice_name]
cice_dg0 = cice[y_lo:y_hi:y_lo, x_lo:x_hi:x_lo, dg0]
print(f"cice_dg0={cice_dg0}")
cice_min = 0
cice_max = 1
if np.any(cice_dg0 <= cice_min) or np.any(cice_dg0 > cice_max):
    print(f"Error: ice concentration at at least one sample point is outside the range 0 ≤ c_ice ≤ 1.")
    sys.exit(8)

# Test the land mask
mask_name = "mask"
mask = data_group[mask_name]
print(f"mask={mask[y_lo:y_hi:y_lo, x_lo:x_hi:x_lo]}")

if mask[50,50] != 1:
    print(f"Error: Test point [50,50] is not wet.")
    sys.exit(8)

if mask[100,140] != 0:
    print(f"Error: Test point [100,140] is not dry.")
    sys.exit(8)

print("Success: integration test passed.")
