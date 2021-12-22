/*!
 * @file DevGrid.cpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGrid.hpp"

#include "/opt/home/include/ncDim.h" // FIXME Remove me
#include "/opt/home/include/ncDouble.h" // FIXME Remove me
#include <ncDim.h>
#include <ncDouble.h>
#include <vector>

// See https://isocpp.org/wiki/faq/pointers-to-members#macro-for-ptr-to-memfn
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

namespace Nextsim {

const std::string DevGrid::ourStructureName = "devgrid";
const std::string xDimName = "x";
const std::string yDimName = "y";
const int DevGrid::nx = 10;
// See https://isocpp.org/wiki/faq/pointers-to-members#array-memfnptrs
// clang-format off
const std::map<std::string, DevGrid::ProgDoubleFn> DevGrid::variableFunctions
    = {    { "hice", &PrognosticData::iceThickness },
           { "cice", &PrognosticData::iceConcentration },
           { "hsnow", &PrognosticData::snowThickness },
           { "sst", &PrognosticData::seaSurfaceTemperature },
           { "sss", &PrognosticData::seaSurfaceSalinity } };
// clang-format on

DevGrid::DevGrid()
    : processedStructureName(ourStructureName)
{
}

DevGrid::~DevGrid() { }

void DevGrid::init(netCDF::NcGroup& metaGroup) { }

void DevGrid::dumpMeta(netCDF::NcGroup& metaGroup) const
{
    // Run the base class output
    IStructure::dumpMeta(metaGroup);
}

void DevGrid::dumpData(netCDF::NcGroup& dataGroup) const
{

    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(xDimName, nx);
    netCDF::NcDim yDim = dataGroup.addDim(yDimName, nx);

    std::vector<netCDF::NcDim> dims = { xDim, yDim };
    for (auto fnNamePair : variableFunctions) {
        std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims));
    }
}

} /* namespace Nextsim */
