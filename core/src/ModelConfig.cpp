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
#include "include/ModelConfigMapElements.ipp"
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
