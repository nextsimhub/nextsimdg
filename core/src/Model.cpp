/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/DevStep.hpp"
#include "include/IDiagnosticOutput.hpp"
#include "include/MissingData.hpp"
#include "include/Module.hpp"
#include "include/StructureFactory.hpp"

#include <string>

// TODO Replace with real logging
#include <iostream>

namespace Nextsim {

const std::string Model::restartOptionName = "model.init_file";

template <>
const std::map<int, std::string> Configured<Model>::keyMap = {
    { Model::RESTARTFILE_KEY, Model::restartOptionName },
#ifdef USE_MPI
    { Model::PARTITIONFILE_KEY, "model.partition_file" },
#endif
    { Model::STARTTIME_KEY, "model.start" },
    { Model::STOPTIME_KEY, "model.stop" },
    { Model::RUNLENGTH_KEY, "model.run_length" },
    { Model::TIMESTEP_KEY, "model.time_step" },
    { Model::MISSINGVALUE_KEY, "model.missing_value" },
};

#ifdef USE_MPI
Model::Model(MPI_Comm comm)
#else
Model::Model()
#endif
{
#ifdef USE_MPI
    m_etadata.setMpiMetadata(comm);
#endif

    iterator.setIterant(&modelStep);

    finalFileName = "restart.nc";
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
    // Configure logging
    Logged::configure();

    startTimeStr = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    durationStr = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    TimePoint timeNow = iterator.parseAndSet(startTimeStr, stopTimeStr, durationStr, stepStr);
    m_etadata.setTime(timeNow);

    MissingData::value
        = Configured::getConfiguration(keyMap.at(MISSINGVALUE_KEY), MissingData::defaultValue);

    initialFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());

    pData.configure();

    modelStep.init();
    modelStep.setInitFile(initialFileName);

#ifdef USE_MPI
    std::string partitionFile
        = Configured::getConfiguration(keyMap.at(PARTITIONFILE_KEY), std::string("partition.nc"));
#endif

#ifdef USE_MPI
    ModelState initialState(
        StructureFactory::stateFromFile(initialFileName, partitionFile, m_etadata));
#else
    ModelState initialState(StructureFactory::stateFromFile(initialFileName));
#endif

    // Get the coordinates from the ModelState for persistence
    m_etadata.extractCoordinates(initialState);

    modelStep.setData(pData);
    modelStep.setMetadata(m_etadata);
    pData.setData(initialState.data);
}

ConfigMap Model::getConfig() const
{
    ConfigMap cMap = {
        { keyMap.at(STARTTIME_KEY), startTimeStr },
        { keyMap.at(STOPTIME_KEY), stopTimeStr },
        { keyMap.at(RUNLENGTH_KEY), durationStr },
        { keyMap.at(TIMESTEP_KEY), stepStr },
        { keyMap.at(MISSINGVALUE_KEY), MissingData::value },
    };

    return cMap;
}

Model::HelpMap& Model::getHelpText(HelpMap& map, bool getAll)
{
    map["Model"] = {
        { keyMap.at(STARTTIME_KEY), ConfigType::STRING, {}, "", "",
            "Model start time, formatted as an ISO8601 date. "
            "Non-calendretical runs can start from year 0 or 1. " },
        { keyMap.at(STOPTIME_KEY), ConfigType::STRING, {}, "", "",
            "Model stop time, formatted as an ISO8601 data. "
            " Will be overridden if a model run length is set. " },
        { keyMap.at(RUNLENGTH_KEY), ConfigType::STRING, {}, "", "",
            "Model run length, formatted as an ISO8601 duration (P prefix). "
            "Overrides the stop time if set. " },
        { keyMap.at(TIMESTEP_KEY), ConfigType::STRING, {}, "", "",
            "Model physics timestep, formatted a ISO8601 duration (P prefix). " },
        { keyMap.at(RESTARTFILE_KEY), ConfigType::STRING, {}, "", "",
            "The file path to the restart file to use for the run." },
        { keyMap.at(MISSINGVALUE_KEY), ConfigType::NUMERIC, { "-∞", "∞" }, "-2³⁰⁰", "",
            "Missing data indicator used for input and output." },
#ifdef USE_MPI
        { keyMap.at(PARTITIONFILE_KEY), ConfigType::STRING, {}, "", "",
            "The file path to the file describing MPI domain decomposition to use for the run." },
#endif
    };

    return map;
}

Model::HelpMap& Model::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    PrognosticData::getHelpRecursive(map, getAll);
    Module::getHelpRecursive<IDiagnosticOutput>(map, getAll);
    return map;
}

void Model::run() { iterator.run(); }

void Model::writeRestartFile()
{
    // TODO Replace with real logging
    Logged::notice(std::string("  Writing state-based restart file: ") + finalFileName + '\n');
    // Copy the configuration from the ModelState to the ModelMetadata
    ConfigMap modelConfig = getConfig();
    modelConfig.merge(pData.getStateRecursive(true).config);
    modelConfig.merge(ConfiguredModule::getAllModuleConfigurations());
    m_etadata.setConfig(modelConfig);
    // Get the model state from PrognosticData and add the coordinates.
    ModelState state = pData.getState();
    m_etadata.affixCoordinates(state);
    StructureFactory::fileFromState(state, m_etadata, finalFileName, true);
}

ModelMetadata& Model::metadata() { return m_etadata; }
} /* namespace Nextsim */
