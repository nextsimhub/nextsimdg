/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/CommonRestartMetadata.hpp"

#include <cstdint>
#include <cstring>
#include <ncInt64.h>
#include <ncString.h>
#include <ncVar.h>

namespace Nextsim {

netCDF::NcFile& CommonRestartMetadata::writeStructureType(
    netCDF::NcFile& ncFile, const ModelMetadata& metadata)
{
    ncFile.putAtt(IStructure::structureNodeName(), metadata.structureName());
    return ncFile;
}

netCDF::NcFile& CommonRestartMetadata::writeRestartMetadata(
    netCDF::NcFile& ncFile, const ModelMetadata& metadata)
{
    // Structure type
    ncFile.putAtt(IStructure::structureNodeName(), metadata.structureName());

    // As Unix time
    netCDF::NcVar unixVar = ncFile.addVar(unformattedName(), netCDF::ncInt64);
    Duration sinceEpoch = metadata.time() - TimePoint();
    std::uint64_t secondsSinceEpoch = sinceEpoch.seconds();
    unixVar.putVar(&secondsSinceEpoch);
    unixVar.putAtt(std::string("units"), "seconds since 1970-01-01T00:00:00Z");

    // Add formatted string as attribute as NetCDF4 does not support string variables
    // in parallel mode
    unixVar.putAtt(std::string("format"), TimePoint::ymdhmsFormat);
    unixVar.putAtt(formattedName(), metadata.m_time.format());

    for (auto entry : metadata.m_config) {
        switch (entry.second.index()) {
        case (ConfigMapType::DOUBLE): {
            ncFile.putAtt(entry.first, netCDF::ncDouble, *std::get_if<double>(&entry.second));
            break;
        }
        case (ConfigMapType::UNSIGNED): {
            ncFile.putAtt(entry.first, netCDF::ncUint, *std::get_if<unsigned>(&entry.second));
            break;
        }
        case (ConfigMapType::INT): {
            ncFile.putAtt(entry.first, netCDF::ncInt, *std::get_if<int>(&entry.second));
            break;
        }
        case (ConfigMapType::STRING): {
            std::string extring = std::get<std::string>(entry.second);
            ncFile.putAtt(entry.first, std::get<std::string>(entry.second));
            break;
        }
        }
    }

    return ncFile;
}

} /* namespace Nextsim */
