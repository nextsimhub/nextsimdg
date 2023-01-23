/*!
 * @file SimpleOutput.cpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SimpleOutput.hpp"

#include "include/Logged.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/StructureFactory.hpp"

#include <sstream>

namespace Nextsim {

void SimpleOutput::outputState(const ModelMetadata& meta)
{
    std::stringstream startStream;
    startStream << meta.time();
    std::string timeFileName = m_filePrefix + "." + startStream.str() + ".nc";
    Logged::info("Outputting "
        + std::to_string(protectedArrayNames.size() + sharedArrayNames.size()) + " fields to "
        + timeFileName + "\n");

    // Create the output by iterating over all fields referenced in ModelState
    ModelState state;
    for (const auto& entry : protectedArrayNames) {
        ModelArrayConstReference macr = getProtectedArray().at(static_cast<size_t>(entry.second));
        if (macr) state.data[entry.first] = *macr;
    }
    for (const auto& entry : sharedArrayNames) {
        ModelArrayReference mar = getSharedArray().at(static_cast<size_t>(entry.second));
        if (mar) state.data[entry.first] = *mar;
    }
    StructureFactory::fileFromState(state, meta, timeFileName);
}
} /* namespace Nextsim */
