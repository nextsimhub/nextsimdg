/*!
 * @file CommonRestartMetadata.cpp
 *
 * @brief The source file for the CommonRestartMetadata class.
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
    // Add a group to the root at the name defined by IStructure.
    netCDF::NcGroup structGroup = rootGroup.addGroup(IStructure::structureNodeName());
    // Add an attribute containing the typeNodeName and the structure name in the metadata.
    structGroup.putAtt(IStructure::typeNodeName(), metadata.structureName());
    return rootGroup;
}

netCDF::NcGroup& CommonRestartMetadata::writeRestartMetadata(
    netCDF::NcGroup& metaGroup, const ModelMetadata& metadata)
{
    // Write the structure type in the metadata group.
    metaGroup.putAtt(IStructure::typeNodeName(), metadata.structureName());

    // Current time is stored in its own sub group, with two entries:
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

    // Loop over all other configuration data and write it to the configuration
    // node within the metadata node.
    netCDF::NcGroup configGroup = metaGroup.addGroup(configurationNode());

    for (auto entry : metadata.m_config) {
        // The netCDF library doesn't have a useful concept of an "any" type,
        // so switch to the correct method to write the attribute into the data
        // file structure.
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
