/*!
 * @file ModelConfig.hpp
 *
 * @date 14 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELCONFIG_HPP
#define MODELCONFIG_HPP

#include "include/ModelMetadata.hpp"
#include "include/Time.hpp"

#include <map>
#include <string>

namespace Nextsim {

class Model;

class ModelConfig {
public:
    const static std::map<int, std::string> keyMap;
    static ConfigMap getConfig();

    // Configuration option that holds the restart file name
    const static std::string restartOptionName;

    enum {
        DUMMY_KEY,
        STARTTIME_KEY,
        STOPTIME_KEY,
        RUNLENGTH_KEY,
        TIMESTEP_KEY,
        RESTARTPERIOD_KEY,
    };

private:
    // Cached values of the start-step-stop/duration times
    static std::string startTimeStr;
    static std::string stopTimeStr;
    static std::string durationStr;
    static std::string stepStr;

    // Period between restart file outputs
    static Duration restartPeriod;

    friend Model;
};

} /* namespace Nextsim */

#endif /* MODELCONFIG_HPP */
