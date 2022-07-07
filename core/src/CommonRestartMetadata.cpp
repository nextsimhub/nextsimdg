/*!
 * @file CommonRestartMetadata.cpp
 *
 * @date Jun 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/CommonRestartMetadata.hpp"

#include <cstdint>
#include <cstring>
#include <ncInt64.h>
#include <ncString.h>
#include <ncVar.h>

#include <iostream> // FIXME remove me

namespace Nextsim {

netCDF::NcGroup& CommonRestartMetadata::writeStructureType(
    netCDF::NcFile& rootGroup, const ModelMetadata& metadata)
{
    netCDF::NcGroup structGroup = rootGroup.addGroup(IStructure::structureNodeName());
    structGroup.putAtt(IStructure::typeNodeName(), metadata.structureName());
    return rootGroup;
}

netCDF::NcGroup& CommonRestartMetadata::writeRestartMetadata(
    netCDF::NcGroup& metaGroup, const ModelMetadata& metadata)
{
    // Structure type
    metaGroup.putAtt(IStructure::typeNodeName(), metadata.structureName());
    // Current time
    netCDF::NcGroup timeGroup = metaGroup.addGroup(timeNodeName());
    // As a formatted string
    netCDF::NcVar formVar = timeGroup.addVar(formattedName(), netCDF::ncString);
    const std::string fTime = metadata.time().format();
    char* timeCopy = new char[fTime.length() + 1];
    std::strcpy(timeCopy, fTime.c_str());
    formVar.putVar(&timeCopy);
    formVar.putAtt(std::string("format"), TimePoint::ymdhmsFormat);
    // As Unix time
    netCDF::NcVar unixVar = timeGroup.addVar(unformattedName(), netCDF::ncInt64);
    Duration sinceEpoch = metadata.time() - TimePoint();
    std::uint64_t secondsSinceEpoch = sinceEpoch.seconds();
    unixVar.putVar(&secondsSinceEpoch);
    unixVar.putAtt(std::string("units"), "seconds since 1970-01-01T00:00:00Z");
    return metaGroup;
}

} /* namespace Nextsim */
