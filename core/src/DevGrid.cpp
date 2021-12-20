/*!
 * @file DevGrid.cpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGrid.hpp"

namespace Nextsim {

const std::string DevGrid::ourStructureName = "devgrid";

DevGrid::DevGrid()
    : processedStructureName(ourStructureName)
{ }

DevGrid::~DevGrid()
{
}

void DevGrid::init(netCDF::NcGroup& metaGroup)
{

}

void DevGrid::dumpMeta(netCDF::NcGroup& metaGroup) const
{
    // Run the base class output
    IStructure::dumpMeta(metaGroup);

}

void DevGrid::dumpData(netCDF::NcGroup& dataGroup) const
{

}


} /* namespace Nextsim */
