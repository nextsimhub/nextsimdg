/*!
 * @file SimpleOutput.cpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SimpleOutput.hpp"
#include "include/Logged.hpp"
#include "include/StructureFactory.hpp"

#include <iostream>
#include <sstream>

namespace Nextsim {

void SimpleOutput::outputState(const ModelState& state, const ModelMetadata& meta)
{
    std::stringstream startStream;
    startStream << meta.time();
    std::string timeFileName = m_filePrefix + "." + startStream.str() + ".nc";
    // Some MPI-IO implemenetations does not like colon in file names
    std::replace(timeFileName.begin(), timeFileName.end(), ':', '_');
    Logged::info("Outputting " + std::to_string(state.data.size()) + " fields and "
        + std::to_string(state.config.size()) + " configurations to " + timeFileName + "\n");
    //    std::cout << "Outputting " << state.size() << " fields to " << timeFileName << std::endl;

    // Copy the configuration from the ModelState to the ModelMetadata
    ModelMetadata metaPlusConfig(meta);
    metaPlusConfig.setConfig(state.config);
    // Create the output
    StructureFactory::fileFromState(state, metaPlusConfig, timeFileName);
}
} /* namespace Nextsim */
