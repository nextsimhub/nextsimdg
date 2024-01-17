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

/*!
 * A class to store model configuration that is to be written to the restart file.
 */
class ModelConfig {
public:
    const static std::map<int, std::string> keyMap;
    static ConfigMap getConfig();

    enum {
        DUMMY_KEY,
        STARTTIME_KEY,
        STOPTIME_KEY,
        RUNLENGTH_KEY,
        TIMESTEP_KEY,
        MISSINGVALUE_KEY,
    };

private:
    // Cached values of the start-step-stop/duration times
    static std::string startTimeStr;
    static std::string stopTimeStr;
    static std::string durationStr;
    static std::string stepStr;

    friend Model;
};

} /* namespace Nextsim */

#endif /* MODELCONFIG_HPP */
