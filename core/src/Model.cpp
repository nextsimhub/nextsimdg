/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/SimpleIterant.hpp"

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
    iterant = new SimpleIterant();
    iterator.setIterant(iterant);

    const int runLength = 5;

    Iterator::TimePoint now(std::chrono::system_clock::now());
    Iterator::Duration dt = std::chrono::seconds(1);
    Iterator::TimePoint hence = now + runLength * dt;

    iterator.setStartStopStep(now, hence, dt);


}

// TODO: add another constructor which takes arguments specifying the
// environment and configuration. This will be the location of the
// logic with selects the components of the model that will run,
// translates I/O details from file configuration to object variable
// values, specifies file paths and likely more besides.

Model::~Model()
{
    if (iterant)
        delete iterant;
}

void Model::configure()
{
    std::string startTimeStr = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    std::string stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    std::string durationStr = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    std::string stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    iterator.parseAndSet(startTimeStr, stopTimeStr, durationStr, stepStr);
}

void Model::run() { iterator.run(); }
} /* namespace Nextsim */
