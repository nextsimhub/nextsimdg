#!/usr/bin/env python3
"""
Convert DG fields in a NetCDF file to the reordered DG basis.

The ordering of the second order DG basis functions was changed from

    1, X, Y, X^2-1/12, Y^2-1/12, XY, ...        (old)

to

    1, X, Y, XY, X^2-1/12, Y^2-1/12, ...        (new)

with X = x-1/2, Y = y-1/2 on the reference element. The components 0-2 and
6-8 are unchanged. Only fields with at least 6 components along the
``dg_comp`` or ``dgstress_comp`` dimension are affected (DG2 advection with 6
components, CG2 stresses with 8 components). Files from the DG1 dynamics are
copied unchanged.

Usage::

    python3 convert_dg_basis_order.py old_restart.nc new_restart.nc
    python3 convert_dg_basis_order.py --reverse new_restart.nc old_restart.nc
"""

import argparse
import sys

import netCDF4
import numpy as np

DG_DIMENSIONS = ("dg_comp", "dgstress_comp")

# new component k is old component OLD_TO_NEW[k]
OLD_TO_NEW = [0, 1, 2, 5, 3, 4, 6, 7, 8]
# old component k is new component NEW_TO_OLD[k]
NEW_TO_OLD = [0, 1, 2, 4, 5, 3, 6, 7, 8]


def permutation(n_comp, reverse):
    """
    Return the index permutation along a DG dimension of size n_comp.

    :param n_comp: number of DG components
    :param reverse: convert from the new to the old ordering
    :return: the permuted indices, or None if no permutation is needed
    """
    if n_comp < 6:
        return None
    perm = NEW_TO_OLD if reverse else OLD_TO_NEW
    return perm[:n_comp]


def copy_group(src, dst, reverse, path=""):
    """
    Recursively copy a NetCDF group, permuting the DG components.

    :param src: source group
    :param dst: destination group
    :param reverse: convert from the new to the old ordering
    :param path: path of the group, for the log
    """
    dst.setncatts(src.__dict__)
    for name, dim in src.dimensions.items():
        dst.createDimension(name, None if dim.isunlimited() else len(dim))

    for name, var in src.variables.items():
        fill = var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None
        out = dst.createVariable(name, var.datatype, var.dimensions, fill_value=fill)
        out.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != "_FillValue"})
        data = var[...]
        for axis, dim_name in enumerate(var.dimensions):
            if dim_name in DG_DIMENSIONS:
                perm = permutation(var.shape[axis], reverse)
                if perm is not None:
                    data = np.take(data, perm, axis=axis)
                    print(f"permuted {path}/{name} along {dim_name}")
        out[...] = data

    for name, grp in src.groups.items():
        copy_group(grp, dst.createGroup(name), reverse, f"{path}/{name}")


def main():
    """Parse the arguments and convert the file."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("input", help="NetCDF file in the old DG basis ordering")
    parser.add_argument("output", help="NetCDF file to write in the new ordering")
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="convert from the new to the old ordering instead",
    )
    args = parser.parse_args()

    if args.input == args.output:
        sys.exit("input and output must be different files")

    with (
        netCDF4.Dataset(args.input, "r") as src,
        netCDF4.Dataset(args.output, "w", format=src.data_model) as dst,
    ):
        src.set_auto_mask(False)
        copy_group(src, dst, args.reverse)


if __name__ == "__main__":
    main()
