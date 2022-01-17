/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/DevGrid.hpp"
#include "include/DevStep.hpp"
#include <string>

namespace Nextsim {

template<>
const std::map<int, std::string> Configured<Model>::keyMap = {
        { Model::RESTARTFILE_KEY, "model.init_file" },
        { Model::STARTTIME_KEY, "model.start" },
        { Model::STOPTIME_KEY, "model.stop" },
        { Model::RUNLENGTH_KEY, "model.run_length" },
        { Model::TIMESTEP_KEY, "model.time_step" },
};

Model::Model()
{
    iterator.setIterant(&iterant);

    dataStructure = nullptr;
}

Model::~Model()
{
}

void Model::configure()
{
    std::string startTimeStr = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    std::string stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    std::string durationStr = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    std::string stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    iterator.parseAndSet(startTimeStr, stopTimeStr, durationStr, stepStr);

    std::string restartFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());

    dataStructure = new DevGrid;
}

void Model::run() { iterator.run(); }
} /* namespace Nextsim */
