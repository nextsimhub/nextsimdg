/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/DevGrid.hpp"
#include "include/DevStep.hpp"
#include "include/DummyExternalData.hpp"
#include "include/StructureFactory.hpp"

#include <string>

// TODO Replace with real logging
#include <iostream>

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

    dataStructure = nullptr;

    finalFileName = "output.nc";
}

Model::~Model()
{
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
    std::string startTimeStr
        = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    std::string stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    std::string durationStr = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    std::string stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    iterator.parseAndSet(startTimeStr, stopTimeStr, durationStr, stepStr);

    initialFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());

    modelStep.setInitFile(initialFileName);

    // Currently, initialize the data here in Model and pass the pointer to the
    // data structure to IModelStep
    dataStructure = StructureFactory::generateFromFile(initialFileName);
    dataStructure->init(initialFileName);
    modelStep.setInitialData(*dataStructure);

    // TODO Real external data handling (in the model step?)
    DummyExternalData::setAll(*dataStructure);
}

void Model::run() { iterator.run(); }

void Model::writeRestartFile()
{
    if (dataStructure) {
        // TODO Replace with real logging
        std::cout << "  Writing restart file: " << finalFileName << std::endl;
        dataStructure->dump(finalFileName);
    }
}
} /* namespace Nextsim */
