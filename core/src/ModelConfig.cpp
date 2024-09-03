/*!
 * @file ModelConfig.cpp
 *
 * @date 14 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelConfig.hpp"
#include "include/MissingData.hpp"

namespace Nextsim {

std::string ModelConfig::startTimeStr;
std::string ModelConfig::stopTimeStr;
std::string ModelConfig::durationStr;
std::string ModelConfig::stepStr;

const std::map<int, std::string> ModelConfig::keyMap = {
    { ModelConfig::STARTTIME_KEY, "model.start" },
    { ModelConfig::STOPTIME_KEY, "model.stop" },
    { ModelConfig::RUNLENGTH_KEY, "model.run_length" },
    { ModelConfig::TIMESTEP_KEY, "model.time_step" },
    { ModelConfig::MISSINGVALUE_KEY, "model.missing_value" },
};

ConfigMap ModelConfig::getConfig()
{
    ConfigMap cMap = {
        { keyMap.at(STARTTIME_KEY), startTimeStr },
        { keyMap.at(STOPTIME_KEY), stopTimeStr },
        { keyMap.at(RUNLENGTH_KEY), durationStr },
        { keyMap.at(TIMESTEP_KEY), stepStr },
        { keyMap.at(MISSINGVALUE_KEY), MissingData::value() },
    };
    return cMap;
}
} /* namespace Nextsim */
