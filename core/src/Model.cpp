/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/DevGrid.hpp"
#include "include/DevStep.hpp"
#include "include/IDiagnosticOutput.hpp"
#include "include/MissingData.hpp"
#include "include/StructureFactory.hpp"

#include <string>

// TODO Replace with real logging
#include <iostream>

namespace Nextsim {

const std::string Model::restartOptionName = "model.init_file";

template <>
const std::map<int, std::string> Configured<Model>::keyMap = {
    { Model::RESTARTFILE_KEY, Model::restartOptionName },
    { Model::PARTITIONFILE_KEY, "model.part_file" },
    { Model::STARTTIME_KEY, "model.start" },
    { Model::STOPTIME_KEY, "model.stop" },
    { Model::RUNLENGTH_KEY, "model.run_length" },
    { Model::TIMESTEP_KEY, "model.time_step" },
};

#ifdef USE_MPI
Model::Model(MPI_Comm comm)
{
    mpiComm = comm;
    CHECK_MPI(MPI_Comm_size(comm, &mpiSize));
    CHECK_MPI(MPI_Comm_rank(comm, &mpiRank));
    iterator.setIterant(&modelStep);
    finalFileName = "restart_" + std::to_string(mpiRank) + ".nc";
}
#else
Model::Model()
{
    iterator.setIterant(&modelStep);

    finalFileName = "restart.nc";
}
#endif // USE_MPI

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

    MissingData mdi;
    mdi.configure();

    initialFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());
    pData.configure();

    modelStep.init();
    modelStep.setInitFile(initialFileName);

#ifdef USE_MPI
    // Read partitioning data from file and load the correct initial state for this process
    std::string partitionFileName
        = Configured::getConfiguration(keyMap.at(PARTITIONFILE_KEY), std::string());
    ModelState initialState(
        StructureFactory::stateFromFile(initialFileName, partitionFileName, mpiComm));
#else
    ModelState initialState(StructureFactory::stateFromFile(initialFileName));
#endif // USE_MPI
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
    };
    // MissingData has a static getState
    cMap.merge(MissingData::getConfig());

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
    };

    return map;
}

Model::HelpMap& Model::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    MissingData::getHelpRecursive(map, getAll);
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

    StructureFactory::fileFromState(pData.getState(), m_etadata, finalFileName);
}

ModelMetadata& Model::metadata() { return m_etadata; }
} /* namespace Nextsim */
