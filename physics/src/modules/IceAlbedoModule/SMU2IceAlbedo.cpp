/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/SMU2IceAlbedo.hpp"
#include "include/KernelAlternatives.hpp"

#include <cmath>

namespace Nextsim {

static constexpr double I0_DEFAULT = 0.17;

double SMU2IceAlbedo::i0;

static const std::string pfx = "SMU2IceAlbedo";
static const std::string i0Key = pfx + ".i0";

std::string SMU2IceAlbedo::getName() const {
    return pfx;
}

void SMU2IceAlbedo::configure() { i0 = Configured::getConfiguration(i0Key, I0_DEFAULT); }

ConfigMap SMU2IceAlbedo::getConfiguration() const
{
    return {
        { i0Key, i0 },
    };
}

/* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
 * alb_ice = 0.64 and alb_sn = 0.85 */

static constexpr double ICE_ALBEDO = 0.64;
static constexpr double SNOW_ALBEDO = 0.85;

void SMU2IceAlbedo::update(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& iceAlbedo = iceAlbedoAccessor.getAutoRW(execSpace);
    auto& icePenSW = icePenSWAccessor.getAutoRW(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);
    const auto& hsnow = hsnowAccessor.getAutoRO(execSpace);

    const double i0 = SMU2IceAlbedo::i0;

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        const double snowThickness = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;

        if (snowThickness > 0.) {
            iceAlbedo[i] = Utils::fmin(
                SNOW_ALBEDO, ICE_ALBEDO + (SNOW_ALBEDO - ICE_ALBEDO) * snowThickness / 0.4);
            icePenSW[i] = 0;
        } else {
            iceAlbedo[i] = ICE_ALBEDO;
            icePenSW[i] = i0;
        }
    });
}
}
