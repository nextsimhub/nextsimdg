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
    const std::string fTime = metadata.m_time.format();
    const char* timeCopy = fTime.c_str();
    formVar.putVar(&timeCopy);
    formVar.putAtt(std::string("format"), TimePoint::ymdhmsFormat);
    // As Unix time
    netCDF::NcVar unixVar = timeGroup.addVar(unformattedName(), netCDF::ncInt64);
    Duration sinceEpoch = metadata.time() - TimePoint();
    std::uint64_t secondsSinceEpoch = sinceEpoch.seconds();
    unixVar.putVar(&secondsSinceEpoch);
    unixVar.putAtt(std::string("units"), "seconds since 1970-01-01T00:00:00Z");

    // All other configuration data
    netCDF::NcGroup configGroup = metaGroup.addGroup(configurationNode());

    for (auto entry : metadata.m_config) {
        switch (entry.second.index()) {
        case (CONFIGMAP_DOUBLE): {
            netCDF::NcVar dblVar = configGroup.addVar(entry.first, netCDF::ncDouble);
            dblVar.putVar(std::get_if<double>(&entry.second));
            break;
        }
        case (CONFIGMAP_UNSIGNED): {
            netCDF::NcVar uintVar = configGroup.addVar(entry.first, netCDF::ncUint);
            uintVar.putVar(std::get_if<unsigned>(&entry.second));
            break;
        }
        case (CONFIGMAP_INT): {
            netCDF::NcVar intVar = configGroup.addVar(entry.first, netCDF::ncInt);
            intVar.putVar(std::get_if<int>(&entry.second));
            break;
        }
        case (CONFIGMAP_STRING): {
            netCDF::NcVar strVar = configGroup.addVar(entry.first, netCDF::ncString);
            std::string extring = std::get<std::string>(entry.second);
            const char* ctring = extring.c_str();
            strVar.putVar(&ctring);
            break;
        }
        }
    }

    return metaGroup;
}

} /* namespace Nextsim */
