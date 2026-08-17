set debuginfod enabled off


# This will add virtual environment if detected
python
import os, glob, site

venv = os.environ.get("VIRTUAL_ENV")
if venv:
    # Locate site-packages inside the active virtual environment
    site_packages = glob.glob(os.path.join(venv, "lib", "python3.*", "site-packages"))
    if site_packages:
        site.addsitedir(site_packages[0])
        print(f"[gdbinit] Loaded virtual environment: {venv}")
end

python

import eigengdb
eigengdb.register_eigen_printers(None)

end


# gdb auto-loading fails due to spack
# Hand-load pretty printers for the C++ standard library
source /usr/share/gcc/python/libstdcxx/v6/printers.py
python register_libstdcxx_printers(None)

source ./gdb/display_marray.py
