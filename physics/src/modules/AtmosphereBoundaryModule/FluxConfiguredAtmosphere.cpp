/*!
 * @file FluxConfiguredAtmosphere.cpp
 *
 * @date Sep 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FluxConfiguredAtmosphere.hpp"

namespace Nextsim {

double FluxConfiguredAtmosphere::qia0 = 305.288;
double FluxConfiguredAtmosphere::dqia_dt0 = 4.5036;
double FluxConfiguredAtmosphere::qow0 = 307.546;
double FluxConfiguredAtmosphere::subl0 = 0;
double FluxConfiguredAtmosphere::snowfall0 = 0;
double FluxConfiguredAtmosphere::rain0 = 0;
double FluxConfiguredAtmosphere::evap0 = 0;
double FluxConfiguredAtmosphere::u0 = 0;
double FluxConfiguredAtmosphere::v0 = 0;

static const std::string pfx = "FluxConfiguredAtmosphere";
static const std::string qiaKey = pfx + ".Q_ia";
static const std::string dqiaKey = pfx + ".dQia_dT";
static const std::string qowKey = pfx + ".Q_ow";
static const std::string sublKey = pfx + ".sublim";
static const std::string snowKey = pfx + ".snowfall";
static const std::string rainKey = pfx + ".rainfall";
static const std::string evapKey = pfx + ".evaporation";
static const std::string uKey = pfx + ".wind_u";
static const std::string vKey = pfx + ".wind_v";

static const std::map<int, std::string> keyMap = {
    { FluxConfiguredAtmosphere::QIA_KEY, qiaKey },
    { FluxConfiguredAtmosphere::DQIA_DT_KEY, dqiaKey },
    { FluxConfiguredAtmosphere::QOW_KEY, qowKey },
    { FluxConfiguredAtmosphere::SUBL_KEY, sublKey },
    { FluxConfiguredAtmosphere::SNOW_KEY, snowKey },
    { FluxConfiguredAtmosphere::RAIN_KEY, rainKey },
    { FluxConfiguredAtmosphere::EVAP_KEY, evapKey },
    { FluxConfiguredAtmosphere::WINDU_KEY, uKey },
    { FluxConfiguredAtmosphere::WINDV_KEY, vKey },
};

ConfigurationHelp::HelpMap& FluxConfiguredAtmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { qiaKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(qia0), "",
            "Total ice to atmosphere heat flux (W m⁻²)." },
        { dqiaKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(dqia_dt0), "",
            "Derivative of the ice atmosphere heat flux with respect to temperature (W m⁻² K⁻¹)." },
        { qowKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(qow0), "",
            "Total open water to atmosphere heat flux (W m⁻²)." },
        { sublKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(subl0), "",
            "Sublimation mass flux from snow to vapour (kg s⁻¹ m⁻²)." },
        { snowKey, ConfigType::NUMERIC, { "0", "∞" }, std::to_string(snowfall0), "",
            "Snowfall mass flux (kg s⁻¹ m⁻²)." },
        { rainKey, ConfigType::NUMERIC, { "0", "∞" }, std::to_string(rain0), "",
            "Rainfall mass flux (kg s⁻¹ m⁻²)." },
        { evapKey, ConfigType::NUMERIC, { "0", "∞" }, std::to_string(subl0), "",
            "Evaporation mass flux (kg s⁻¹ m⁻²)." },
        { uKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(u0), "",
            "Component of wind in the x (eastward) direction (m s⁻¹)." },
        { vKey, ConfigType::NUMERIC, { "-∞", "∞" }, std::to_string(v0), "",
            "Component of wind in the y (northward) direction (m s⁻¹)." },
    };
    return map;
}
void FluxConfiguredAtmosphere::configure()
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

void FluxConfiguredAtmosphere::setData(const ModelState::DataMap& dm)
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
