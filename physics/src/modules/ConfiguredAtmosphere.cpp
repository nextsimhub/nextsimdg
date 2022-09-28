/*!
 * @file ConfiguredAtmosphere.cpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredAtmosphere.hpp"

namespace Nextsim {

double ConfiguredAtmosphere::qia0 = 305.288;
double ConfiguredAtmosphere::dqia_dt0 = 4.5036;
double ConfiguredAtmosphere::qow0 = 307.546;
double ConfiguredAtmosphere::subl0 = 0;
double ConfiguredAtmosphere::snowfall0 = 0;
double ConfiguredAtmosphere::rain0 = 0;
double ConfiguredAtmosphere::evap0 = 0;
double ConfiguredAtmosphere::u0 = 0;
double ConfiguredAtmosphere::v0 = 0;

template <>
const std::map<int, std::string> Configured<ConfiguredAtmosphere>::keyMap = {
    { ConfiguredAtmosphere::QIA_KEY, "ConfiguredAtmosphere.Q_ia" },
    { ConfiguredAtmosphere::DQIA_DT_KEY, "ConfiguredAtmosphere.dQia_dT" },
    { ConfiguredAtmosphere::QOW_KEY, "ConfiguredAtmosphere.Q_ow" },
    { ConfiguredAtmosphere::SUBL_KEY, "ConfiguredAtmosphere.sublim" },
    { ConfiguredAtmosphere::SNOW_KEY, "ConfiguredAtmosphere.snowfall" },
    { ConfiguredAtmosphere::RAIN_KEY, "ConfiguredAtmosphere.rainfall" },
    { ConfiguredAtmosphere::WINDU_KEY, "ConfiguredAtmosphere.wind_u" },
    { ConfiguredAtmosphere::WINDV_KEY, "ConfiguredAtmosphere.wind_v" },
};

void ConfiguredAtmosphere::configure()
{
    qia0 = Configured::getConfiguration(keyMap.at(QIA_KEY), qia0);
    dqia_dt0 = Configured::getConfiguration(keyMap.at(DQIA_DT_KEY), dqia_dt0);
    qow0 = Configured::getConfiguration(keyMap.at(QOW_KEY), qow0);
    subl0 = Configured::getConfiguration(keyMap.at(SUBL_KEY), subl0);
    snowfall0 = Configured::getConfiguration(keyMap.at(SNOW_KEY), snowfall0);
    rain0 = Configured::getConfiguration(keyMap.at(RAIN_KEY), rain0);
    evap0 = Configured::getConfiguration(keyMap.at(EVAP_KEY), evap0);
    u0 = Configured::getConfiguration(keyMap.at(WINDU_KEY), u0);
    v0 = Configured::getConfiguration(keyMap.at(WINDV_KEY), v0);
}

void ConfiguredAtmosphere::setData(const ModelState::DataMap& dm)
{
    IAtmosphereBoundary::setData(dm);
    qia = qia0;
    dqia_dt = dqia_dt0;
    qow = qow0;
    subl = subl0;
    snow = snowfall0;
    rain = rain0;
    evap = evap0;
    uwind = u0;
    vwind = v0;
}

} /* namespace Nextsim */
