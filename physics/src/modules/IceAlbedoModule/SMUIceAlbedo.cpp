/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/SMUIceAlbedo.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

static constexpr FloatType I0_DEFAULT = 0.17;

FloatType SMUIceAlbedo::i0;

static const std::string pfx = "SMUIceAlbedo";
static const std::string i0Key = pfx + ".i0";

std::string SMUIceAlbedo::getName() const { return pfx; }

void SMUIceAlbedo::configure() { i0 = Configured::getConfiguration(i0Key, I0_DEFAULT); }

ConfigMap SMUIceAlbedo::getConfiguration() const
{
    return {
        { i0Key, i0 },
    };
}

/* This scheme mimics Semtner 76 and Maykut and Untersteiner 71 when
 * alb_ice = 0.64 and alb_sn = 0.85 */

static constexpr FloatType ICE_ALBEDO = 0.64;
static constexpr FloatType SNOW_ALBEDO = 0.85;

void SMUIceAlbedo::update(const TimestepTime& tst)
{
    auto execSpace = DefaultExecutionSpace();
    auto& iceAlbedo = iceAlbedoAccessor.getAutoRW(execSpace);
    auto& icePenSW = icePenSWAccessor.getAutoRW(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);
    const auto& hsnow = hsnowAccessor.getAutoRO(execSpace);

    const FloatType i0 = SMUIceAlbedo::i0;

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        const FloatType snowThickness = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;

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
