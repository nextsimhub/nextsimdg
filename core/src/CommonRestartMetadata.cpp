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
    // As Unix time
    netCDF::NcVar unixVar = timeGroup.addVar(unformattedName(), netCDF::ncInt64);
    Duration sinceEpoch = metadata.time() - TimePoint();
    std::uint64_t secondsSinceEpoch = sinceEpoch.seconds();
    unixVar.putVar(&secondsSinceEpoch);
    unixVar.putAtt(std::string("units"), "seconds since 1970-01-01T00:00:00Z");

    // Add formatted string as attribute as NetCDF4 does not support string variables
    // in parallel mode
    unixVar.putAtt(std::string("format"), TimePoint::ymdhmsFormat);
    unixVar.putAtt(formattedName(), metadata.m_time.format());

    // All other configuration data
    netCDF::NcGroup configGroup = metaGroup.addGroup(configurationNode());

    for (auto entry : metadata.m_config) {
        switch (entry.second.index()) {
        case (ConfigMapType::DOUBLE): {
            configGroup.putAtt(entry.first, netCDF::ncDouble, *std::get_if<double>(&entry.second));
            break;
        }
        case (ConfigMapType::UNSIGNED): {
            configGroup.putAtt(entry.first, netCDF::ncUint, *std::get_if<unsigned>(&entry.second));
            break;
        }
        case (ConfigMapType::INT): {
            configGroup.putAtt(entry.first, netCDF::ncInt, *std::get_if<int>(&entry.second));
            break;
        }
        case (ConfigMapType::STRING): {
            std::string extring = std::get<std::string>(entry.second);
            configGroup.putAtt(entry.first, std::get<std::string>(entry.second));
            break;
        }
        }
    }

    return metaGroup;
}

} /* namespace Nextsim */
