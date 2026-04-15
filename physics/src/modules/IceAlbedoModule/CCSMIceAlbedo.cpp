/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/CCSMIceAlbedo.hpp"
#include "include/KernelAlternatives.hpp"

#include <cmath>

/* Albedo scheme from ccsm 3 */
/* The scheme is simplified using the assumption that visible solar
 * radiation accounts for 52% of the spectrum and the near-infrared for
 * 48% (the same split as used in cice when running in stand-alone
 * mode). */

/* This is the ccsm3 scheme when alb_ice = 0.538 and alb_sn = 0.8256 */

namespace Nextsim {

static constexpr double ICE_ALBEDO0 = 0.538;
static constexpr double SNOW_ALBEDO0 = 0.8256;
static constexpr double I0_DEFAULT = 0.17;

double CCSMIceAlbedo::iceAlbedo = ICE_ALBEDO0;
double CCSMIceAlbedo::snowAlbedo = SNOW_ALBEDO0;
double CCSMIceAlbedo::i0 = I0_DEFAULT;

static const std::string pfx = "CCSMIceAlbedo";
static const std::string iceAlbedoKey = pfx + ".iceAlbedo";
static const std::string snowAlbedoKey = pfx + ".snowAlbedo";
static const std::string i0Key = pfx + ".I0";

std::string CCSMIceAlbedo::getName() const { return pfx; }

void CCSMIceAlbedo::configure()
{
    iceAlbedo = Configured::getConfiguration(iceAlbedoKey, ICE_ALBEDO0);
    snowAlbedo = Configured::getConfiguration(snowAlbedoKey, SNOW_ALBEDO0);
    i0 = Configured::getConfiguration(i0Key, I0_DEFAULT);
}

ConfigMap CCSMIceAlbedo::getConfiguration() const
{
    return {
        { iceAlbedoKey, iceAlbedo },
        { snowAlbedoKey, snowAlbedo },
        { i0Key, i0 },
    };
}

void CCSMIceAlbedo::update(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& iceAlbedo = iceAlbedoAccessor.getAutoRW(execSpace);
    auto& icePenSW = icePenSWAccessor.getAutoRW(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);
    const auto& hsnow = hsnowAccessor.getAutoRO(execSpace);
    const auto& tsurf = tsurfAccessor.getAutoRO(execSpace);

    const double iceAlbedoBase = CCSMIceAlbedo::iceAlbedo;
    const double snowAlbedoBase = CCSMIceAlbedo::snowAlbedo;
    const double i0 = CCSMIceAlbedo::i0;

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        const double snowThickness = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;
        const double temperature = tsurf[i];

        constexpr double tLimit = -1.;
        const double iceAlbedoT = iceAlbedoBase - Utils::fmax(0., 0.075 * (temperature - tLimit));
        const double snowAlbedoT = snowAlbedoBase - Utils::fmax(0., 0.124 * (temperature - tLimit));
        const double snowCoverFraction = snowThickness / (snowThickness + 0.02);

        const double albedo
            = snowCoverFraction * snowAlbedoT + (1 - snowCoverFraction) * iceAlbedoT;
        const double penSW = (1. - snowCoverFraction) * i0;

        iceAlbedo[i] = albedo;
        icePenSW[i] = penSW;
    });
}

CCSMIceAlbedo::HelpMap& CCSMIceAlbedo::getHelpText(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { iceAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, ConfigurationHelp::toString(ICE_ALBEDO0),
            "", "Albedo of snow-free ice." },
        { snowAlbedoKey, ConfigType::NUMERIC, { "0", "1" },
            ConfigurationHelp::toString(SNOW_ALBEDO0), "", "Albedo of snow." },
        // From the former documentation of IIceAlbedo::surfaceShortWaveBalance:
        // Fraction of short-wave radiation that can penetrate bare ice (not taking snow cover into
        // account)
        { i0Key, ConfigType::NUMERIC, { "0", "1" }, ConfigurationHelp::toString(I0_DEFAULT), "",
            "Transmissivity of ice." },

    };
    return map;
}
CCSMIceAlbedo::HelpMap& CCSMIceAlbedo::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}
} /* namespace Nextsim */
