/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/DevGrid.hpp"
#include "include/DevStep.hpp"
#include "include/MissingData.hpp"
#include "include/ModelState.hpp"
#include "include/StructureFactory.hpp"

#include <string>

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<Model>::keyMap = {
    { Model::RESTARTFILE_KEY, "model.init_file" },
    { Model::STARTTIME_KEY, "model.start" },
    { Model::STOPTIME_KEY, "model.stop" },
    { Model::RUNLENGTH_KEY, "model.run_length" },
    { Model::TIMESTEP_KEY, "model.time_step" },
};

Model::Model()
{
    iterator.setIterant(&modelStep);

    finalFileName = "restart.nc";
}

Model::~Model()
{
    // Perform any required shutdown of the coupler, if included.
    coupler.terminate();
    /*
     * Try writing out a valid restart file. If the model and computer are in a
     * state where this can be completed, great! If they are not then the
     * restart file is unlikely to be valid or otherwise stored properly, and
     * we abandon the writing.
     */
    try {
        writeRestartFile();
    } catch (std::exception& e) {
        // If there are any exceptions at all, fail without writing
    }
}

void Model::configure()
{
    // Configure logging
    Logged::configure();

    std::string startTimeStr
        = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    std::string stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    std::string durationStr = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    std::string stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    TimePoint timeNow = iterator.parseAndSet(startTimeStr, stopTimeStr, durationStr, stepStr);
    data.metadata.setTime(timeNow);

    MissingData mdi;
    mdi.configure();

    // Configure the coupler, if included.
    coupler.configure();

    initialFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());

    data.pData.configure();

    modelStep.init();
    modelStep.setInitFile(initialFileName);

    ModelState initialState(StructureFactory::stateFromFile(initialFileName));
    data.pData.setData(initialState);
    modelStep.setData(data);

    data.aoState.configure();

    // With everything configured, do the metadata dependent initialization of
    // the coupler, if included.
    coupler.initialize(data.metadata);
}

void Model::run() {
    iterator.run();
}

void Model::writeRestartFile()
{
    // TODO Replace with real logging
    Logged::notice(std::string("  Writing state-based restart file: ") + finalFileName + '\n');
    StructureFactory::fileFromState(data.pData.getState(), data.getMetadata(), finalFileName);
}

} /* namespace Nextsim */
