/*!
 * @file Model.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Model.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/DevGrid.hpp"
#include "include/DevStep.hpp"
#include "include/IDiagnosticOutput.hpp"
#include "include/MissingData.hpp"
#include "include/Module.hpp"
#include "include/StructureFactory.hpp"

#include <string>

// TODO Replace with real logging
#include <iostream>

namespace Nextsim {

template <typename Key, typename Value>
static std::map<Key, Value> combineMaps(std::map<Key, Value>& map1, std::map<Key, Value>& map2)
{
    std::map<Key, Value> map3 = map1;
    map3.merge(map2);
    return map3;
}

static std::map<int, std::string> interimKeyMap = ModelConfig::keyMap;
const std::string Model::restartOptionName = "model.init_file";
static std::map<int, std::string> restartKeyMap = {
    { Model::RESTARTFILE_KEY, Model::restartOptionName },
    { Model::RESTARTPERIOD_KEY, "model.restart_period" },
    { Model::RESTARTOUTFILE_KEY, "model.restart_file" },
};

template <>
const std::map<int, std::string> Configured<Model>::keyMap
    = combineMaps(interimKeyMap, restartKeyMap);

Model::Model()
{
    iterator.setIterant(&modelStep);

    finalFileName = std::string("restart") + TimePoint::ymdhmsFormat + ".nc";
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

    ModelConfig::startTimeStr
        = Configured::getConfiguration(keyMap.at(STARTTIME_KEY), std::string());
    ModelConfig::stopTimeStr = Configured::getConfiguration(keyMap.at(STOPTIME_KEY), std::string());
    ModelConfig::durationStr
        = Configured::getConfiguration(keyMap.at(RUNLENGTH_KEY), std::string());
    ModelConfig::stepStr = Configured::getConfiguration(keyMap.at(TIMESTEP_KEY), std::string());

    TimePoint timeNow = iterator.parseAndSet(ModelConfig::startTimeStr, ModelConfig::stopTimeStr,
        ModelConfig::durationStr, ModelConfig::stepStr);
    m_etadata.setTime(timeNow);

    MissingData::value
        = Configured::getConfiguration(keyMap.at(MISSINGVALUE_KEY), MissingData::defaultValue);

    initialFileName = Configured::getConfiguration(keyMap.at(RESTARTFILE_KEY), std::string());
    finalFileName = Configured::getConfiguration(keyMap.at(RESTARTOUTFILE_KEY), finalFileName);

    pData.configure();

    modelStep.init();
    modelStep.setInitFile(initialFileName);

    ModelState initialState(StructureFactory::stateFromFile(initialFileName));

    std::string restartPeriodStr
        = Configured::getConfiguration(keyMap.at(RESTARTPERIOD_KEY), std::string());
    restartPeriod = Duration(restartPeriodStr);
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
        { keyMap.at(RESTARTPERIOD_KEY), ConfigType::STRING, {}, "", "",
            "The period between restart file outputs, formatted as an ISO8601 duration (P prefix) "
            "or number of seconds." },
        { keyMap.at(MISSINGVALUE_KEY), ConfigType::NUMERIC, { "-∞", "∞" }, "-2³⁰⁰", "",
            "Missing data indicator used for input and output." },
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
    std::string formattedFileName = m_etadata.time().format(finalFileName);

    Logged::notice(std::string("  Writing state-based restart file: ") + formattedFileName + '\n');
    // Copy the configuration from the ModelState to the ModelMetadata
    ConfigMap modelConfig = getConfig();
    modelConfig.merge(pData.getStateRecursive(true).config);
    modelConfig.merge(ConfiguredModule::getAllModuleConfigurations());
    m_etadata.setConfig(modelConfig);
    StructureFactory::fileFromState(pData.getState(), m_etadata, formattedFileName, true);
}

ModelMetadata& Model::metadata() { return m_etadata; }
} /* namespace Nextsim */
