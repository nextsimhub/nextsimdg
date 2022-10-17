/*!
 * @file DGNetcdf.cpp
 *
 * @date Oct 4, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DGNetcdf.hpp"

#include <ncFile.h>

namespace Nextsim {

DGNetcdf::DGNetcdf() { p_file = new netCDF::NcFile(); }

DGNetcdf::DGNetcdf(const std::string& filePath, FileMode fm)
{
    netCDF::NcFile::FileMode modeMux
        = { netCDF::NcFile::FileMode::read, netCDF::NcFile::FileMode::write,
              netCDF::NcFile::FileMode::replace, netCDF::NcFile::FileMode::newFile };
    p_file = new netCDF::NcFile(filePath, modeMux[fm]);
}

DGNetcdf::~DGNetcdf()
{
    if (p_file) {
        p_file->close();
        delete p_file;
    }
}

netCDF::NcGroup& DGNetcdf::write(netCDF::NcGroup& group, const ModelArray& dgField, const std::string& fieldName)
{

    return group;
}


} /* namespace Nextsim */
