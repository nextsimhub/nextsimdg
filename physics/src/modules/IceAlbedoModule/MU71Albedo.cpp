/*!
 * @author  Einar Örn Ólason <einar.olason@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/MU71Albedo.hpp"
#include "include/KernelAlternatives.hpp"
#include <array>

namespace Nextsim {

static constexpr FloatType I0_DEFAULT = 0.17;

FloatType MU71Albedo::i0;

static const std::string pfx = "MU71Albedo";
static const std::string i0Key = pfx + ".i0";

std::string MU71Albedo::getName() const { return pfx; }

void MU71Albedo::configure() { i0 = Configured::getConfiguration(i0Key, I0_DEFAULT); }

ConfigMap MU71Albedo::getConfiguration() const
{
    return {
        { i0Key, i0 },
    };
}

// Monthly snow albedo from Maykut and Untersteiner (1971)
static constexpr std::array<FloatType, 12> ALBEDO_TABLE
    = { 0.85, 0.85, 0.83, 0.81, 0.82, 0.78, 0.64, 0.69, 0.84, 0.85, 0.85, 0.85 };
//      Jan,  Feb,  Mar,  Apr,  Mai,  Jun,  Jul,  Aug,  Sept, Oct,  Nov,  Dec

MU71Albedo::MU71Albedo()
    : snowAlbedo(std::vector(ALBEDO_TABLE.begin(), ALBEDO_TABLE.end()))
{
}

static constexpr FloatType ICE_ALBEDO = 0.64;

void MU71Albedo::update(const TimestepTime& tst)
{
    // always do update on the host because it currently not worth the effort to port
    // monthlyCubicBSpline
    auto& iceAlbedo = iceAlbedoAccessor.getHostRW();
    auto& icePenSW = icePenSWAccessor.getHostRW();
    const auto& cice = ciceAccessor.getHostRO();
    const auto& hsnow = hsnowAccessor.getHostRO();

    const TimePoint& tp = tst.start;
    const FloatType dayOfYear = tp.gmtime()->tm_yday;
    const bool isLeap = ((tp.gmtime()->tm_year % 4 == 0) && (tp.gmtime()->tm_year % 100 != 0))
        || (tp.gmtime()->tm_year % 400 == 0);
    const FloatType i0 = MU71Albedo::i0;

    overElements([&](const size_t i) {
        const FloatType snowThickness = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;

        if (snowThickness == 0.0_ft) {
            iceAlbedo[i] = ICE_ALBEDO;
            icePenSW[i] = i0;
        } else {
            iceAlbedo[i] = snowAlbedo(dayOfYear, isLeap);
            icePenSW[i] = 0.;
        }
    });
}

}
