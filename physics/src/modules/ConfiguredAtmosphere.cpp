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

static const std::string pfx = "ConfiguredAtmosphere";
static const std::string qiaKey = pfx + ".Q_ia";
static const std::string dqiaKey = pfx + ".dQia_dT";
static const std::string qowKey = pfx + ".Q_ow";
static const std::string sublKey = pfx + ".sublim";
static const std::string snowKey = pfx + ".snowfall";
static const std::string rainKey = pfx + ".rainfall";
static const std::string uKey = pfx + ".wind_u";
static const std::string vKey = pfx + ".wind_v";

template <>
const std::map<int, std::string> Configured<ConfiguredAtmosphere>::keyMap = {
    { ConfiguredAtmosphere::QIA_KEY, qiaKey },
    { ConfiguredAtmosphere::DQIA_DT_KEY, dqiaKey },
    { ConfiguredAtmosphere::QOW_KEY, qowKey },
    { ConfiguredAtmosphere::SUBL_KEY, sublKey },
    { ConfiguredAtmosphere::SNOW_KEY, snowKey },
    { ConfiguredAtmosphere::RAIN_KEY, rainKey },
    { ConfiguredAtmosphere::WINDU_KEY, uKey },
    { ConfiguredAtmosphere::WINDV_KEY, vKey },
};

ConfigurationHelp::HelpMap& ConfiguredAtmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
            { qiaKey, ConfigType::NUMERIC, {"-∞", "∞"}, std::to_string(qia0), "", "Total ice to atmosphere heat flux." },
            { dqiaKey, ConfigType::NUMERIC, {"-∞", "∞"}, std::to_string(dqia_dt0), "", "Derivative of the ice atmosphere heat flux with respect to temperature." },
            { qowKey, ConfigType::NUMERIC, {"-∞", "∞"}, std::to_string(qow0), "", "Total open water to atmosphere heat flux." },
            { sublKey, ConfigType::NUMERIC, {"-∞", "∞"}, std::to_string(subl0), "", "Sublimation mass flux from snow to vapour." },
            { snowKey, ConfigType::NUMERIC, {"0", "∞"}, std::to_string(snowfall0), "", "Snowfall mass flux." },
            { rainKey, ConfigType::NUMERIC, {"0", "∞"}, std::to_string(rain0), "", "Rainfall mass flux." },
            { uKey, ConfigType::NUMERIC, {"-∞", "∞"}, std::to_string(u0), "", "Component of wind in the x (eastward) direction." },
            { vKey, ConfigType::NUMERIC, {"-∞", "∞"}, std::to_string(v0), "", "Component of wind in the y (northward) direction." },
    };
    return map;
}
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
