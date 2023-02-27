#include <ncCheck.h>
#include <netcdf.h>
#include <netcdf_par.h>

#include "include/ParallelNetcdfFile.hpp"

using namespace netCDF;

// Not sure why it is needed but let's replicate netCDF::ncFile::open
// in this respect
extern int g_ncid;

NcFilePar::NcFilePar(
    const std::string& filePath, const FileMode fMode, MPI_Comm comm, MPI_Info mpiInfo)
{
    open_par(filePath, fMode, comm, mpiInfo);
}

void NcFilePar::open_par(
    const std::string& filePath, const FileMode fMode, MPI_Comm comm, MPI_Info mpiInfo)
{
    if (!nullObject)
        close();

    switch (fMode) {
    case NcFile::write:
        ncCheck(nc_open_par(filePath.c_str(), NC_WRITE, comm, mpiInfo, &myId), __FILE__, __LINE__);
        break;
    case NcFile::read:
        ncCheck(
            nc_open_par(filePath.c_str(), NC_NOWRITE, comm, mpiInfo, &myId), __FILE__, __LINE__);
        break;
    case NcFile::newFile:
        ncCheck(nc_create_par(filePath.c_str(), NC_NETCDF4 | NC_NOCLOBBER, comm, mpiInfo, &myId),
            __FILE__, __LINE__);
        break;
    case NcFile::replace:
        ncCheck(nc_create_par(filePath.c_str(), NC_NETCDF4 | NC_CLOBBER, comm, mpiInfo, &myId),
            __FILE__, __LINE__);
        break;
    }

    g_ncid = myId;

    nullObject = false;
}
