/*!
 * @file ParaGridIO.cpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParaGridIO.hpp"

#include "include/gridNames.hpp"

#include <ncFile.h>
#include <ncGroup.h>

#include <map>
#include <string>

namespace Nextsim {

ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    ModelState state;
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    // Dimensions
    std::multimap<std::string, netCDF::NcDim> dimMap = dataGroup.getDims();

    for (auto entry : ModelArray::definedDimensions) {
        netCDF::NcDim dim = dataGroup.getDim(entry.second.name);
        if (dim.isNull()) {
            throw std::out_of_range(std::string("No dimension named ") + entry.second.name + " found in " + filePath);
        }
        entry.second.length = dim.getSize();
    }



    return state;
}
} /* namespace Nextsim */
