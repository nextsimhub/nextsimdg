/*!
 * @file FluxConfiguredAtmosphere.cpp
 *
 * @date 30 Apr 2025
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
        { qiaKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(qia0), "",
            "Total ice to atmosphere heat flux (W m⁻²)." },
        { dqiaKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(dqia_dt0), "",
            "Derivative of the ice atmosphere heat flux with respect to temperature (W m⁻² K⁻¹)." },
        { qowKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(qow0), "",
            "Total open water to atmosphere heat flux (W m⁻²)." },
        { sublKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(subl0), "",
            "Sublimation mass flux from snow to vapour (kg s⁻¹ m⁻²)." },
        { snowKey, ConfigType::NUMERIC, { "0", "∞" }, ConfigurationHelp::toString(snowfall0), "",
            "Snowfall mass flux (kg s⁻¹ m⁻²)." },
        { rainKey, ConfigType::NUMERIC, { "0", "∞" }, ConfigurationHelp::toString(rain0), "",
            "Rainfall mass flux (kg s⁻¹ m⁻²)." },
        { evapKey, ConfigType::NUMERIC, { "0", "∞" }, ConfigurationHelp::toString(subl0), "",
            "Evaporation mass flux (kg s⁻¹ m⁻²)." },
        { uKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(u0), "",
            "Component of wind in the x (eastward) direction (m s⁻¹)." },
        { vKey, ConfigType::NUMERIC, { "-∞", "∞" }, ConfigurationHelp::toString(v0), "",
            "Component of wind in the y (northward) direction (m s⁻¹)." },
    };
    return map;
}
void FluxConfiguredAtmosphere::configure()
{
    qia0 = getConfiguration(keyMap.at(QIA_KEY), qia0);
    dqia_dt0 = getConfiguration(keyMap.at(DQIA_DT_KEY), dqia_dt0);
    qow0 = getConfiguration(keyMap.at(QOW_KEY), qow0);
    subl0 = getConfiguration(keyMap.at(SUBL_KEY), subl0);
    snowfall0 = getConfiguration(keyMap.at(SNOW_KEY), snowfall0);
    rain0 = getConfiguration(keyMap.at(RAIN_KEY), rain0);
    evap0 = getConfiguration(keyMap.at(EVAP_KEY), evap0);
    u0 = getConfiguration(keyMap.at(WINDU_KEY), u0);
    v0 = getConfiguration(keyMap.at(WINDV_KEY), v0);
}

ConfigMap FluxConfiguredAtmosphere::getConfiguration() const
{
    return {
       { keyMap.at(QIA_KEY), qia0 },
       { keyMap.at(DQIA_DT_KEY), dqia_dt0 },
       { keyMap.at(QOW_KEY), qow0 },
       { keyMap.at(SUBL_KEY), subl0 },
       { keyMap.at(SNOW_KEY), snowfall0 },
       { keyMap.at(RAIN_KEY), rain0 },
       { keyMap.at(EVAP_KEY), evap0 },
       { keyMap.at(WINDU_KEY), u0 },
       { keyMap.at(WINDV_KEY), v0 },
    };
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

void FluxConfiguredAtmosphere::update(const TimestepTime& tst)
{
    /* The open water heat flux is reset by the thermodynamics, so that the (slab) ocean doesn't
     * super cool. Therefore, we have to reset it here on every update */
    qow = qow0;
}

} /* namespace Nextsim */
