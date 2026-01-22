/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SMUIceAlbedo.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

double SMUIceAlbedo::m_i0;

static constexpr double i0Default = 0.17;

static const std::unordered_map<int, std::string> keyMap = {
    { SMUIceAlbedo::I0_KEY, "nextsim_thermo.I_0" },
};

void SMUIceAlbedo::configure()
{
    m_i0 = Configured::getConfiguration(keyMap.at(I0_KEY), i0Default);
}

ConfigMap SMUIceAlbedo::getConfiguration() const
{
    return {
        { keyMap.at(I0_KEY), m_i0 },
    };
}

/* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
 * alb_ice = 0.64 and alb_sn = 0.85 */

const double ICE_ALBEDO = 0.64;
const double SNOW_ALBEDO = 0.85;

std::tuple<double, double> SMUIceAlbedo::surfaceShortWaveBalance(
    double temperature, double snowThickness, double i0)
{
    double albedo, penSW;
    if (snowThickness > 0.) {
        albedo = SNOW_ALBEDO;
        penSW = 0.;
    } else {
        albedo = ICE_ALBEDO;
        penSW = i0;
    }
    return { albedo, penSW };
}

void SMUIceAlbedo::update(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& iceAlbedo = iceAlbedoAccessor.getAutoRW(execSpace);
    auto& icePenSW = icePenSWAccessor.getAutoRW(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);
    const auto& hsnow = hsnowAccessor.getAutoRO(execSpace);

    const double i0 = SMUIceAlbedo::m_i0;

    overElements(OVER_ELEMENTS_LAMBDA(ElementIndex i) {
        const double snowThickness = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;

        if (snowThickness > 0.) {
            iceAlbedo[i] = SNOW_ALBEDO;
            icePenSW[i] = 0.;
        } else {
            iceAlbedo[i] = ICE_ALBEDO;
            icePenSW[i] = i0;
        }
    });
}

}
