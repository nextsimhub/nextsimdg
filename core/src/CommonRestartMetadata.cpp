/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/CommonRestartMetadata.hpp"
#include "include/ModelMetadata.hpp"

#include "include/MissingData.hpp"

#include <cstdint>
#include <cstring>
#include <ncDouble.h>
#include <ncFloat.h>
#include <ncInt64.h>
#include <ncString.h>
#include <ncVar.h>

namespace Nextsim {

netCDF::NcFile& CommonRestartMetadata::writeStructureType(netCDF::NcFile& ncFile)
{
    auto& metadata = ModelMetadata::getInstance();
    ncFile.putAtt(IStructure::structureNodeName(), metadata.structureName());
    return ncFile;
}

netCDF::NcFile& CommonRestartMetadata::writeRestartMetadata(netCDF::NcFile& ncFile)
{
    // Structure type
    auto& metadata = ModelMetadata::getInstance();
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

    ncFile.putAtt(startTimeName(), metadata.startTime().format());
    ncFile.putAtt(stopTimeName(), metadata.stopTime().format());
    ncFile.putAtt(stepLengthName(), std::to_string(metadata.stepLength().seconds()));
    ncFile.putAtt(runLengthName(), std::to_string(metadata.runLength().seconds()));

    ncFile.putAtt(missingDataName(), netCDF::NcDouble(), MissingData::value());

    for (auto entry : metadata.m_config) {
        switch (entry.second.index()) {
        case (ConfigMapType::DOUBLE): {
            ncFile.putAtt(entry.first, netCDF::ncDouble, *std::get_if<double>(&entry.second));
            break;
        }
        case (ConfigMapType::FLOAT): {
            ncFile.putAtt(entry.first, netCDF::ncFloat, *std::get_if<float>(&entry.second));
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
