/*!
 * @file MinimumIce.cpp
 *
 * @date Nov 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/MinimumIce.hpp"

namespace Nextsim {

double MinimumIce::minc;
double MinimumIce::minh;

static const double mincDefault = 1e-12;
static const double minhDefault = 0.01;

template <>
const std::map<int, std::string> Configured<MinimumIce>::keyMap = {
    { MinimumIce::MINC_KEY, "nextsim_thermo.min_conc" },
    { MinimumIce::MINH_KEY, "nextsim_thermo.min_thick" },

};

MinimumIce::HelpMap& MinimumIce::getHelpText(HelpMap& map, bool getAll)
{
    map["MinimumIce"] = {
        { keyMap.at(MINC_KEY), ConfigType::NUMERIC, { "0", "1" }, std::to_string(mincDefault), "",
            "Minimum allowed ice concentration." },
        { keyMap.at(MINH_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(minhDefault), "m",
            "Minimum allowed ice thickness." },
    };
    return map;
}
MinimumIce::HelpMap& MinimumIce::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    return map;
}

void MinimumIce::configure()
{
    // Configure constants
    minc = Configured::getConfiguration(keyMap.at(MINC_KEY), mincDefault);
    minh = Configured::getConfiguration(keyMap.at(MINH_KEY), minhDefault);
}

ConfigMap MinimumIce::getConfiguration() const
{
    return {
        { keyMap.at(MINC_KEY), minc },
        { keyMap.at(MINH_KEY), minh },
    };
}

} /* namespace Nextsim */
