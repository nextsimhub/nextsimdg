/*!
 * @file ConfiguredAtmosphere.cpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfiguredAtmosphere.hpp"

#include "include/Module.hpp"

namespace Nextsim {

double ConfiguredAtmosphere::tair0 = 0; // Freezing
double ConfiguredAtmosphere::tdew0 = 0; // Fog
double ConfiguredAtmosphere::pair0 = 101325; // Mean pressure
double ConfiguredAtmosphere::sw0 = 0; // Night
double ConfiguredAtmosphere::lw0 = 315.637; // Stefan-Boltzmann at 0˚C
double ConfiguredAtmosphere::snowfall0 = 0; // No snow
double ConfiguredAtmosphere::rain0 = 0; // No rain
double ConfiguredAtmosphere::windspeed0 = 0; // Still

static const std::string pfx = "ConfiguredAtmosphere";
static const std::string tKey = pfx + ".t_air";
static const std::string tdewKey = pfx + ".t_dew";
static const std::string pKey = pfx + ".pmsl";
static const std::string swKey = pfx + ".sw_in";
static const std::string lwKey = pfx + ".lw_in";
static const std::string snowKey = pfx + ".snowfall";
static const std::string rainKey = pfx + ".rainfall";
static const std::string windKey = pfx + ".wind_speed";

template <>
const std::map<int, std::string> Configured<ConfiguredAtmosphere>::keyMap = {
    { ConfiguredAtmosphere::TAIR_KEY, tKey },
    { ConfiguredAtmosphere::TDEW_KEY, tdewKey },
    { ConfiguredAtmosphere::PAIR_KEY, pKey },
    { ConfiguredAtmosphere::SW_KEY, swKey },
    { ConfiguredAtmosphere::LW_KEY, lwKey },
    { ConfiguredAtmosphere::SNOW_KEY, snowKey },
    { ConfiguredAtmosphere::RAIN_KEY, rainKey },
    { ConfiguredAtmosphere::WIND_KEY, windKey },
};

ConfigurationHelp::HelpMap& ConfiguredAtmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { tKey, ConfigType::NUMERIC, { "-210", "374" }, std::to_string(tair0), "",
            "Air temperature at 2 m (˚C)." },
        { tdewKey, ConfigType::NUMERIC, { "-273", "374" }, std::to_string(tdew0), "",
            "Dew point temperature at 2 m (˚C)." },
        { pKey, ConfigType::NUMERIC, { "0", "3395800" }, std::to_string(pair0), "",
            "Surface air pressure (Pa)." },
        { swKey, ConfigType::NUMERIC, { "0", "1.390e122" }, std::to_string(sw0), "",
            "Downward solar radiation (W m⁻²)." },
        { lwKey, ConfigType::NUMERIC, { "0", "1.390e122" }, std::to_string(lw0), "",
            "Downward non-solar radiation (W m⁻²)." },
        { snowKey, ConfigType::NUMERIC, { "0", "2.998e8" }, std::to_string(snowfall0), "",
            "Snowfall mass flux (kg s⁻¹ m⁻²)." },
        { rainKey, ConfigType::NUMERIC, { "0", "2.998e8" }, std::to_string(rain0), "",
            "Rainfall mass flux (kg s⁻¹ m⁻²)." },
        { windKey, ConfigType::NUMERIC, { "0", "2.998e8" }, std::to_string(windspeed0), "",
            "Windspeed (m s⁻¹)." },
    };
    Module::getHelpRecursive<IFluxCalculation>(map, getAll);

    return map;
}

void ConfiguredAtmosphere::configure()
{
    snowfall0 = Configured::getConfiguration(keyMap.at(SNOW_KEY), snowfall0);
    rain0 = Configured::getConfiguration(keyMap.at(RAIN_KEY), rain0);

    fluxImpl = &Module::getImplementation<IFluxCalculation>();
    tryConfigure(fluxImpl);
}

void ConfiguredAtmosphere::setData(const ModelState::DataMap& dm)
{
    IAtmosphereBoundary::setData(dm);
    snow = snowfall0;
    rain = rain0;
}

void ConfiguredAtmosphere::update(const TimestepTime& tst)
{
    fluxImpl->update(tst);
}

} /* namespace Nextsim */
