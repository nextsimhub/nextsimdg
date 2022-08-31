/*!
 * @file ConfiguredAtmosphere.cpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredAtmosphere.hpp"

namespace Nextsim {

double ConfiguredAtmosphere::tair0 = -1;
double ConfiguredAtmosphere::tdew0 = -0.5;
double ConfiguredAtmosphere::pair0 = 1e5;
double ConfiguredAtmosphere::rmix0 = -1;
double ConfiguredAtmosphere::sw_in0 = 0;
double ConfiguredAtmosphere::lw_in0 = 311;
double ConfiguredAtmosphere::snowfall0 = 0;
double ConfiguredAtmosphere::windspeed0 = 0;

template <>
const std::map<int, std::string> Configured<ConfiguredAtmosphere>::keyMap = {
    { ConfiguredAtmosphere::TAIR_KEY, "ConfiguredAtmosphere.T_air" },
    { ConfiguredAtmosphere::TDEW_KEY, "ConfiguredAtmosphere.T_dew" },
    { ConfiguredAtmosphere::PAIR_KEY, "ConfiguredAtmosphere.p_air" },
    { ConfiguredAtmosphere::RMIX_KEY, "ConfiguredAtmosphere.r_mix" },
    { ConfiguredAtmosphere::SWIN_KEY, "ConfiguredAtmosphere.sw_in" },
    { ConfiguredAtmosphere::LWIN_KEY, "ConfiguredAtmosphere.lw_in" },
    { ConfiguredAtmosphere::SNOW_KEY, "ConfiguredAtmosphere.snowfall" },
    { ConfiguredAtmosphere::WIND_KEY, "ConfiguredAtmosphere.wind_speed" },
};

void ConfiguredAtmosphere::configure()
{
    tair = Configured::getConfiguration(keyMap.at(TAIR_KEY), tair0);
    tdew = Configured::getConfiguration(keyMap.at(TDEW_KEY), tdew0);
    pair = Configured::getConfiguration(keyMap.at(PAIR_KEY), pair0);
    rmix = Configured::getConfiguration(keyMap.at(RMIX_KEY), rmix0);
    sw_in = Configured::getConfiguration(keyMap.at(SWIN_KEY), sw_in0);
    lw_in = Configured::getConfiguration(keyMap.at(LWIN_KEY), lw_in0);
    snowfall = Configured::getConfiguration(keyMap.at(SNOW_KEY), snowfall0);
    windSpeed = Configured::getConfiguration(keyMap.at(WIND_KEY), windspeed0);
}

} /* namespace Nextsim */
