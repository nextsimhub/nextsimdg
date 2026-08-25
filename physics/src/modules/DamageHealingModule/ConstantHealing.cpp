/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/ConstantHealing.hpp"
#include "include/IceMinima.hpp"
#include "include/KernelAlternatives.hpp"
#include "kokkos/include/KokkosTimer.hpp"

namespace Nextsim {

FloatType ConstantHealing::tD = 0.;
static const FloatType tDDefault = 15;

static const std::map<int, std::string> keyMap
    = { { ConstantHealing::TD_KEY, "ConstantHealing.td" } };

void ConstantHealing::configure()
{
    // the option is defined in days, but the model wants seconds
    tD = Configured::getConfiguration(keyMap.at(TD_KEY), tDDefault);
    tD *= 86400.;
}

ConfigMap ConstantHealing::getConfiguration() const { return { { keyMap.at(TD_KEY), tD } }; }

ConstantHealing::HelpMap& ConstantHealing::getHelpText(HelpMap& map, bool getAll)
{
    map["ConstantHealing"] = { { keyMap.at(TD_KEY), ConfigType::NUMERIC, { "0", "∞" },
        ConfigurationHelp::toString(tDDefault), "days",
        "The healing time scale (t_d) for brittle rheologies" } };
    return map;
}

ConstantHealing::HelpMap& ConstantHealing::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

/* Heal damage through
 * 1. Lateral ice formation: Newly formed ice is undamaged
 * 2. Constant healing with a given time scale (tD) */
void ConstantHealing::update(const TimestepTime& tstep)
{
    auto execSpace = DefaultExecutionSpace();
    auto& damage = damageAccessor.getAutoRW(execSpace);
    const auto& deltaCi = deltaCiAccessor.getAutoRO(execSpace);
    const auto& cice = ciceAccessor.getAutoRO(execSpace);

    const FloatType tD = ConstantHealing::tD;
    const FloatType cMin = IceMinima::c();
    const FloatType dt = tstep.step.seconds();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // No ice, no healing
        if (cice[i] <= cMin) {
            damage[i] = 1.;
            return;
        }

        // Only lateral growth contributes to healing, not melt(!)
        FloatType const lateralGrowth = Utils::max(0.0_ft, deltaCi[i]);

        /* 1. Lateral ice formation
         * A weighted average of the original damage, weighted by the old concentration, and the
         * undamaged new ice damage (1), weighted by the concentration of new ice. */
        damage[i] = (damage[i] * (cice[i] - lateralGrowth) + lateralGrowth) / cice[i];

        /* 2. Constant healing
         * Damage healing using a constant timescale. Originally conceived as an exponential
         * decay, but then revised to a linear one. */
        // This is what Sylvain and Pierre (Bouillon and Rampal, 2015)
        // damage[i] +=  damage[i] * tstep.step / tD;

        // This is what Véro did (Dansereau et al., 2016)
        damage[i] += dt / tD;
        damage[i] = Utils::min(1.0_ft, damage[i]);
    });
}

}
