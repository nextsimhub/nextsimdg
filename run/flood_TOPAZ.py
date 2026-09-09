import glob
import os
import sys
from pathlib import Path

import xarray as xr
from flooding import FloodInterpolator

"""
Create a flooded version of TOPAZ input files. Assume we have a directory structure like this:

daily_mean
├── 2006
│   ├── topaz_rean_200610.nc
│   ├── topaz_rean_200611.nc
│   └── topaz_rean_200612.nc
└── 2007
    ├── topaz_rean_200701.nc
    ├── topaz_rean_200702.nc
    ├── topaz_rean_200703.nc
    ├── topaz_rean_200704.nc
    └── topaz_rean_200705.nc

and create this:

daily_mean_flooded
├── 2006
│   ├── topaz_rean_200610.nc
│   ├── topaz_rean_200611.nc
│   └── topaz_rean_200612.nc
└── 2007
    ├── topaz_rean_200701.nc
    ├── topaz_rean_200702.nc
    ├── topaz_rean_200703.nc
    ├── topaz_rean_200704.nc
    └── topaz_rean_200705.nc
"""

if len(sys.argv) != 2 or sys.argv[1] == "--help":
    print("Usage: python " + sys.argv[0] + " <root directory of TOPAZ files>")
    print(
        "  We expect a directory structure where each year is a sub-directory of the root directory given as an argument."
    )
    sys.exit(0)

root_dir = sys.argv[1]

# Create an output directory. Fail if it exists
out_dir = root_dir.rstrip("/") + "_flooded"
os.makedirs(out_dir, exist_ok=False)

# Just pick a random netCDF in "root_dir" to build the FloodInterpolator
ds = xr.load_dataset(next(iter(Path(root_dir).rglob("*.nc"))))
flooder = FloodInterpolator(ds.x.values, ds.y.values, ds.mask.values.astype(bool))
ds.close()

for directory in glob.glob(root_dir + "/*"):
    basedir = os.path.basename(directory)
    os.mkdir(out_dir + "/" + basedir)
    for file in glob.glob(root_dir + "/" + basedir + "/*.nc"):
        ds = xr.open_dataset(file)
        flooder.flood_dataset(ds, ["thetao", "so", "mlotst", "siconc", "sithick", "sisnthick"])
        flooder.flood_zeros(ds, ["vxo", "vyo"])
        flooder.flood_nearest(ds, ["zos"])

        encoding = {
            var: {"zlib": True, "complevel": 4, "shuffle": False} for var in ds.data_vars
        }
        print(os.path.join(out_dir, basedir, os.path.basename(file)))
        ds.to_netcdf(
            os.path.join(out_dir, basedir, os.path.basename(file)), encoding=encoding
        )
        ds.close()
