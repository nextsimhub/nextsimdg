/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/WintonAlbedo.hpp"
#include "include/KernelAlternatives.hpp"

/* Albedo scheme from Winton (2000)
 * It's a very simplified scheme: "snow was given an albedo of 0.8---reduced to 0.75 under melting
 * conditions---and ice an albedo of 0.65. The latter value was tuned to give an annual mean ice
 * thickness of 3 m." */

namespace Nextsim {

static constexpr FloatType ICE_ALBEDO0 = 0.65;
static constexpr FloatType SNOW_ALBEDO0 = 0.80;
static constexpr FloatType MELT_ALBEDO0 = 0.75;
static constexpr FloatType I0_DEFAULT = 0.17;

FloatType WintonAlbedo::iceAlbedo = ICE_ALBEDO0;
FloatType WintonAlbedo::snowAlbedo = SNOW_ALBEDO0;
FloatType WintonAlbedo::meltAlbedo = MELT_ALBEDO0;
FloatType WintonAlbedo::i0 = I0_DEFAULT;

static const std::string pfx = "WintonAlbedo";
static const std::string iceAlbedoKey = pfx + ".iceAlbedo";
static const std::string snowAlbedoKey = pfx + ".snowAlbedo";
static const std::string meltAlbedoKey = pfx + ".meltAlbedo";
static const std::string i0Key = pfx + ".i0";

std::string WintonAlbedo::getName() const { return pfx; }

void WintonAlbedo::configure()
{
    iceAlbedo = Configured::getConfiguration(iceAlbedoKey, ICE_ALBEDO0);
    snowAlbedo = Configured::getConfiguration(snowAlbedoKey, SNOW_ALBEDO0);
    meltAlbedo = Configured::getConfiguration(meltAlbedoKey, MELT_ALBEDO0);
    i0 = Configured::getConfiguration(i0Key, I0_DEFAULT);
}

ConfigMap WintonAlbedo::getConfiguration() const
{
    return {
        { iceAlbedoKey, iceAlbedo },
        { snowAlbedoKey, snowAlbedo },
        { meltAlbedoKey, meltAlbedo },
        { i0Key, i0 },
    };
}

WintonAlbedo::HelpMap& WintonAlbedo::getHelpText(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { iceAlbedoKey, ConfigType::NUMERIC, { "0", "1" }, ConfigurationHelp::toString(ICE_ALBEDO0),
            "", "Albedo of snow-free ice." },
        { snowAlbedoKey, ConfigType::NUMERIC, { "0", "1" },
            ConfigurationHelp::toString(SNOW_ALBEDO0), "", "Albedo of dry snow." },
        { meltAlbedoKey, ConfigType::NUMERIC, { "0", "1" },
            ConfigurationHelp::toString(MELT_ALBEDO0), "", "Albedo of melting snow." },
    };
    return map;
}
WintonAlbedo::HelpMap& WintonAlbedo::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

void WintonAlbedo::update(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& iceAlbedo = iceAlbedoAccessor.getAutoRW(execSpace);
    auto& icePenSW = icePenSWAccessor.getAutoRW(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);
    const auto& hsnow = hsnowAccessor.getAutoRO(execSpace);
    const auto& tsurf = tsurfAccessor.getAutoRO(execSpace);

    const FloatType iceAlbedoBase = WintonAlbedo::iceAlbedo;
    const FloatType snowAlbedoBase = WintonAlbedo::snowAlbedo;
    const FloatType meltAlbedoBase = WintonAlbedo::meltAlbedo;
    const FloatType i0 = WintonAlbedo::i0;

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        const FloatType snowThickness = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;

        if (snowThickness > 0.0_ft) {
            iceAlbedo[i] = tsurf[i] < 0. ? snowAlbedoBase : meltAlbedoBase;
            icePenSW[i] = 0.;
        } else {
            iceAlbedo[i] = iceAlbedoBase;
            icePenSW[i] = i0;
        }
    });
}

} /* namespace Nextsim */
