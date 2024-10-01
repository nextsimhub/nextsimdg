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
#include "include/Finalizer.hpp"
#include "include/IDiagnosticOutput.hpp"
#include "include/MissingData.hpp"
#include "include/NextsimModule.hpp"
#include "include/StructureFactory.hpp"

#include <string>

// TODO Replace with real logging
#include <iostream>

namespace Nextsim {

// Map of configuration that will be written to the restart file
const std::string Model::restartOptionName = "model.init_file";
// Map of all configuration keys for the main model, including those not to be
// written to the restart file.
static const std::map<int, std::string> keyMap = {
#include "include/ModelConfigMapElements.ipp"
    { Model::RESTARTFILE_KEY, Model::restartOptionName },
#ifdef USE_MPI
    { Model::PARTITIONFILE_KEY, "model.partition_file" },
#endif
    { Model::STARTTIME_KEY, "model.start" },
    { Model::STOPTIME_KEY, "model.stop" },
    { Model::RUNLENGTH_KEY, "model.run_length" },
    { Model::TIMESTEP_KEY, "model.time_step" },
    { Model::MISSINGVALUE_KEY, "model.missing_value" },
    { Model::RESTARTPERIOD_KEY, "model.restart_period" },
    { Model::RESTARTOUTFILE_KEY, "model.restart_file" },
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

    finalFileName = std::string("restart") + TimePoint::ymdhmsFormat + ".nc";
}

Model::~Model() { }

void Model::configure()
{
    // Configure logging
    Logged::configure();

    // Store the start/stop/step configuration directly in ModelConfig before
    // parsing these values to the numerical time values used by the model.
    ModelConfig::startTimeStr
        = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    ModelConfig::stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    ModelConfig::durationStr
        = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    ModelConfig::stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    // Set the time correspond to the current (initial) model state
    TimePoint timeNow = iterator.parseAndSet(ModelConfig::startTimeStr, ModelConfig::stopTimeStr,
        ModelConfig::durationStr, ModelConfig::stepStr);
    m_etadata.setTime(timeNow);

    // Configure the missing data value
    MissingData::setValue(
        Configured::getConfiguration(keyMap.at(MISSINGVALUE_KEY), MissingData::defaultValue));

    // Parse the initial restart file name and the pattern for output restart files
    initialFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());
    finalFileName = Configured::getConfiguration(keyMap.at(RESTARTOUTFILE_KEY), finalFileName);

    pData.configure();

    modelStep.init();
    modelStep.setInitFile(initialFileName);

#ifdef USE_MPI
    std::string partitionFile
        = Configured::getConfiguration(keyMap.at(PARTITIONFILE_KEY), std::string("partition.nc"));
    m_etadata.getPartitionMetadata(partitionFile);
#endif

#ifdef USE_MPI
    ModelState initialState(StructureFactory::stateFromFile(initialFileName, m_etadata));
#else
    ModelState initialState(StructureFactory::stateFromFile(initialFileName));
#endif

    // The period with which to write restart files.
    std::string restartPeriodStr
        = Configured::getConfiguration(keyMap.at(RESTARTPERIOD_KEY), std::string("0"));
    restartPeriod = Duration(restartPeriodStr);

    // Get the coordinates from the ModelState for persistence
    m_etadata.extractCoordinates(initialState);

    modelStep.setData(pData);
    modelStep.setMetadata(m_etadata);
    modelStep.setRestartDetails(restartPeriod, finalFileName);
    pData.setData(initialState.data);
}

ConfigMap Model::getConfig() const
{
    ConfigMap cMap = ModelConfig::getConfig();
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
            "Model physics timestep, formatted as an ISO8601 duration (P prefix). " },
        { keyMap.at(RESTARTFILE_KEY), ConfigType::STRING, {}, "", "",
            "The file path to the restart file to use for the run." },
        { keyMap.at(RESTARTPERIOD_KEY), ConfigType::STRING, {}, "0", "",
            "The period between restart file outputs, formatted as an ISO8601 "
            "duration (P prefix) or number of seconds. A value of zero "
            "ensures no intermediate restart files are written." },
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

void Model::run()
{
    iterator.run();
    writeRestartFile();
    Finalizer::finalize();
}

//! Write a restart file for the model.
void Model::writeRestartFile()
{
    std::string formattedFileName = m_etadata.time().format(finalFileName);
    pData.writeRestartFile(formattedFileName, m_etadata);
}

ModelMetadata& Model::metadata() { return m_etadata; }
} /* namespace Nextsim */
