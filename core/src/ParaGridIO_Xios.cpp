/*!
 * @file ParaGridIO_Xios.cpp
 *
 * @date May 1, 2024
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 */

#include "include/ParaGridIO_Xios.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/FileCallbackCloser.hpp"
#include "include/MissingData.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncException.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <algorithm>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>

namespace Nextsim {

ParaGridIO::~ParaGridIO() = default;

ModelState getModelState(const std::string& filePath, ModelMetadata& metadata) override
{
    ModelState state;
    // TODO: Implement this method?
    return state;
}

ModelState ParaGridIO::readForcingTimeStatic(
    const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath)
{
    ModelState state;
    // TODO: Implement this method?
    return state;
}

void ParaGridIO::dumpModelState(
    const ModelState& state, const ModelMetadata& metadata, const std::string& filePath)
{
    // TODO: Implement this method?
}

void ParaGridIO::writeDiagnosticTime(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
{
    // TODO: Implement this method?
}

void ParaGridIO::close(const std::string& filePath)
{
    // TODO: Implement this method?
}

} /* namespace Nextsim */
