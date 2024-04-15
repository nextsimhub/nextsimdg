/*!
 * @file WintonAlbedo.cpp
 *
 * @date Jan 11, 2024
 * @author Tim Spain
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/WintonAlbedo.hpp"

/* Albedo scheme from Winton (2000)
 * It's a very simplified scheme: "snow was given an albedo of 0.8---reduced to 0.75 under melting
 * conditions---and ice an albedo of 0.65. The latter value was tuned to give an annual mean ice
 * thickness of 3 m." */

namespace Nextsim {

static const double ICE_ALBEDO0 = 0.65;
static const double SNOW_ALBEDO0 = 0.80;
static const double MELT_ALBEDO0 = 0.75;

double WintonAlbedo::iceAlbedo = ICE_ALBEDO0;
double WintonAlbedo::snowAlbedo = SNOW_ALBEDO0;
double WintonAlbedo::meltAlbedo = MELT_ALBEDO0;

static const std::string pfx = "WintonAlbedo";
static const std::string iceAlbedoKey = pfx + ".iceAlbedo";
static const std::string snowAlbedoKey = pfx + ".snowAlbedo";
static const std::string meltAlbedoKey = pfx + ".meltAlbedo";

std::tuple<double, double> WintonAlbedo::surfaceShortWaveBalance(
    double temperature, double snowThickness, double i0)
{
    double albedo, penSW;
    if (snowThickness > 0.) {
        penSW = 0.;
        if (temperature < 0.)
            albedo = snowAlbedo;
        else
            albedo = meltAlbedo;
    } else {
        penSW = i0;
        albedo = iceAlbedo;
    }

    return { albedo, penSW };
}

void WintonAlbedo::configure()
{
    iceAlbedo = Configured::getConfiguration(iceAlbedoKey, ICE_ALBEDO0);
    snowAlbedo = Configured::getConfiguration(snowAlbedoKey, SNOW_ALBEDO0);
    meltAlbedo = Configured::getConfiguration(meltAlbedoKey, MELT_ALBEDO0);
}

ConfigMap WintonAlbedo::getConfiguration() const
{
    return {
        { iceAlbedoKey, iceAlbedo },
        { snowAlbedoKey, snowAlbedo },
        { meltAlbedoKey, meltAlbedo },
    };
}

WintonAlbedo::HelpMap& WintonAlbedo::getHelpText(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { iceAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, std::to_string(ICE_ALBEDO0), "",
            "Albedo of snow-free ice." },
        { snowAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, std::to_string(SNOW_ALBEDO0), "",
            "Albedo of dry snow." },
        { meltAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, std::to_string(MELT_ALBEDO0), "",
            "Albedo of melting snow." },
    };
    return map;
}
WintonAlbedo::HelpMap& WintonAlbedo::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}
} /* namespace Nextsim */
