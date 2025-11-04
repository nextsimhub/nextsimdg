/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Kacper Kornet <kk562@cam.ac.uk>
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
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif

#include <string>

// TODO Replace with real logging
#include <iostream>

namespace Nextsim {

// Map of configuration that will be written to the restart file
const std::string Model::restartOptionName = "model.init_file";
// Map of all configuration keys for the main model, including those not to be
// written to the restart file.
static const std::map<int, std::string> keyMap = {
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

Model::Model()
    : iterator(modelStep)
{
#ifdef USE_MPI
    std::string partitionFile
        = Configured::getConfiguration(keyMap.at(PARTITIONFILE_KEY), std::string("partition.nc"));
    auto& metadata = ModelMetadata::getInstance(partitionFile);
#else
    auto& metadata = ModelMetadata::getInstance();
#endif
    metadata.setFinalFileName("restart" + TimePoint::ymdhmsFormat + ".nc");
}

Model::~Model() { }

void Model::configureRestarts()
{
    ModelMetadata& metadata = ModelMetadata::getInstance();

    // Parse the initial restart file name and the pattern for output restart files
    metadata.setInitialFileName(
        Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string()));
    metadata.setFinalFileName(
        Configured::getConfiguration(keyMap.at(RESTARTOUTFILE_KEY), metadata.finalFileName()));

    // The period with which to write restart files.
    std::string restartPeriodStr
        = Configured::getConfiguration(keyMap.at(RESTARTPERIOD_KEY), std::string("0"));
    metadata.setRestartPeriod(Duration(restartPeriodStr));
}

void Model::configureTime()
{
    ModelMetadata& metadata = ModelMetadata::getInstance();

    // Start/stop times. Run length will override stop time, if present.
    std::string startTimeStr
        = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    std::string stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    std::string runLengthStr
        = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    std::string stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    if (runLengthStr.empty()) {
        if (stopTimeStr.empty()) {
            throw std::invalid_argument(std::string("At least one of ") + keyMap.at(STOPTIME_KEY)
                + " or " + keyMap.at(RUNLENGTH_KEY) + " must be set");
        } else {
            metadata.setTimes(startTimeStr, TimePoint(stopTimeStr), stepStr);
        }
    } else {
        metadata.setTimes(startTimeStr, Duration(runLengthStr), stepStr);
    }
    iterator.setStartStopStep(metadata.startTime(), metadata.stopTime(), metadata.stepLength());
}

void Model::configure()
{
    // Configure logging
    Logged::configure();

    // Configure the missing data value
    MissingData::setValue(
        Configured::getConfiguration(keyMap.at(MISSINGVALUE_KEY), MissingData::defaultValue));

    configureRestarts();

    pData.configure();

    auto& metadata = ModelMetadata::getInstance();
    modelStep.init();
    modelStep.setInitFile(metadata.initialFileName());

    configureTime();

    ModelState initialState(StructureFactory::stateFromFile(metadata.initialFileName()));

    // Get the coordinates from the ModelState for persistence
    metadata.extractCoordinates(initialState);

    modelStep.setData(pData);
    modelStep.setRestartDetails(metadata.restartPeriod(), metadata.finalFileName());
    pData.setData(initialState.data);
}

Model::HelpMap& Model::getHelpText(HelpMap& map, bool getAll)
{
    map["Model"] = {
        { keyMap.at(STARTTIME_KEY), ConfigType::STRING, {}, "", "",
            "Model start time, formatted as an ISO8601 date. "
            "Non-calendretical runs can start from year 0 or 1. " },
        { keyMap.at(STOPTIME_KEY), ConfigType::STRING, {}, "", "",
            "Model stop time, formatted as an ISO8601 date. "
            " Will be overridden if a model run length is set. " },
        { keyMap.at(RUNLENGTH_KEY), ConfigType::STRING, {}, "", "",
            "Model run length, formatted as an ISO8601 duration (P prefix) or an integer number of "
            "seconds. "
            "Overrides the stop time if set. " },
        { keyMap.at(TIMESTEP_KEY), ConfigType::STRING, {}, "", "",
            "Model physics timestep, formatted as an ISO8601 duration (P prefix) or an integer "
            "number of seconds. " },
        { keyMap.at(RESTARTFILE_KEY), ConfigType::STRING, {}, "", "",
            "The file path to the initial restart file to be used for the run." },
        { keyMap.at(RESTARTPERIOD_KEY), ConfigType::STRING, {}, "0", "",
            "The period between restart file outputs, formatted as an ISO8601 "
            "duration (P prefix) or number of seconds. A value of zero "
            "ensures no intermediate restart files are written." },
        { keyMap.at(RESTARTOUTFILE_KEY), ConfigType::STRING, {}, "restart%Y-%m-%dT%H:%M:%SZ.nc", "",
            "The file name or file name pattern to be used for the output restart files. A fixed "
            "name can be used, which will be overwritten every time a restart checkpoint is "
            "created. Otherwise the file name may include formatting codes supported by the C/C++ "
            "time library to create a new file at every restart checkpoint." },
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
#ifdef USE_XIOS
    Xios::getHelpRecursive(map, getAll);
#endif
    return map;
}

void Model::run()
{
    try {
        iterator.run();
    } catch (const std::exception& e) {
        writeRestartFile();
        Finalizer::finalize();
        throw;
    }

    writeRestartFile();
    Finalizer::finalize();
}

//! Write a restart file for the model.
void Model::writeRestartFile()
{
    auto& metadata = ModelMetadata::getInstance();
    std::string formattedFileName = metadata.time().format(metadata.finalFileName());
    pData.writeRestartFile(formattedFileName);
}

} /* namespace Nextsim */
