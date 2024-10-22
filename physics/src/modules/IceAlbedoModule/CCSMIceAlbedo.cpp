/*!
 * @file CCSMIceAlbedo.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain
 */

#include "include/CCSMIceAlbedo.hpp"

#include <cmath>

/* Albedo scheme from ccsm 3 */
/* The scheme is simplified using the assumption that visible solar
 * radiation accounts for 52% of the spectrum and the near-infrared for
 * 48% (the same split as used in cice when running in stand-alone
 * mode). */

/* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */

namespace Nextsim {

static const double ICE_ALBEDO0 = 0.538;
static const double SNOW_ALBEDO0 = 0.8256;

double CCSMIceAlbedo::iceAlbedo = ICE_ALBEDO0;
double CCSMIceAlbedo::snowAlbedo = SNOW_ALBEDO0;

static const std::string pfx = "CCSMIceAlbedo";
static const std::string iceAlbedoKey = pfx + ".iceAlbedo";
static const std::string snowAlbedoKey = pfx + ".snowAlbedo";

std::tuple<double, double> CCSMIceAlbedo::surfaceShortWaveBalance(
    double temperature, double snowThickness, double i0)
{
    const double tLimit = -1.;
    double iceAlbedoT = iceAlbedo - std::fmax(0., 0.075 * (temperature - tLimit));
    double snowAlbedoT = snowAlbedo - std::fmax(0., 0.124 * (temperature - tLimit));
    double snowCoverFraction = snowThickness / (snowThickness + 0.02);

    const double albedo = snowCoverFraction * snowAlbedoT + (1 - snowCoverFraction) * iceAlbedoT;
    const double penSW = (1. - snowCoverFraction) * i0;

    return { albedo, penSW };
}

void CCSMIceAlbedo::configure()
{
    iceAlbedo = Configured::getConfiguration(iceAlbedoKey, ICE_ALBEDO0);
    snowAlbedo = Configured::getConfiguration(snowAlbedoKey, SNOW_ALBEDO0);
}

ConfigMap CCSMIceAlbedo::getConfiguration() const
{
    return {
        { iceAlbedoKey, iceAlbedo },
        { snowAlbedoKey, snowAlbedo },
    };
}

CCSMIceAlbedo::HelpMap& CCSMIceAlbedo::getHelpText(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { iceAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, std::to_string(ICE_ALBEDO0), "",
            "Albedo of snow-free ice." },
        { snowAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, std::to_string(SNOW_ALBEDO0), "",
            "Albedo of snow." },
    };
    return map;
}
CCSMIceAlbedo::HelpMap& CCSMIceAlbedo::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}
} /* namespace Nextsim */
