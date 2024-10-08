/*!
 * @file SimpleOutput.cpp
 *
 * @date 24 Sep 2024
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
    // Some MPI-IO implemenetations does not like colon in file names
    std::replace(timeFileName.begin(), timeFileName.end(), ':', '_');
    Logged::info(
        "Outputting " + std::to_string(externalNames.size()) + " fields to " + timeFileName + "\n");

    // Create the output by iterating over all fields referenced in ModelState
    ModelState state;
    auto storeData = ModelComponent::getStore().getAllData();
    for (auto entry : storeData) {
        if (entry.second)
            state.data[entry.first] = *entry.second;
    }
    StructureFactory::fileFromState(state, meta, timeFileName);
}
} /* namespace Nextsim */
