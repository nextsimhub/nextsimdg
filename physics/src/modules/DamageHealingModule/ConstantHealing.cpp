/*!
 * @file ConstantHealing.hpp
 *
 * @date Jun 3, 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/ConstantHealing.hpp"

namespace Nextsim {

double ConstantHealing::tD = 0.;
static const double tDDefault = 15;

static const std::map<int, std::string> localKeyMap
    = { { ConstantHealing::TD_KEY, "ConstantHealing.td" } };

void ConstantHealing::configure()
{
    // the option is defined in days, but the model wants seconds
    tD = Configured::getConfiguration(localKeyMap.at(TD_KEY), tDDefault);
    tD *= 86400.;
}

ModelState ConstantHealing::getStateRecursive(const Nextsim::OutputSpec& os) const
{
    return { {}, { { localKeyMap.at(TD_KEY), tD } } };
}

ConstantHealing::HelpMap& ConstantHealing::getHelpText(HelpMap& map, bool getAll)
{
    map["ConstantHealing"]
        = { { localKeyMap.at(TD_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(tDDefault),
            "days", "The healing time scale (t_d) for brittle rheologies" } };
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
    overElements(std::bind(&ConstantHealing::updateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tstep);
}

void ConstantHealing::updateElement(size_t i, const TimestepTime& tstep)
{
    // Only lateral growth contributes to healing, not melt(!)
    double const lateralGrowth = std::max(0., deltaCi[i]);

    /* 1. Lateral ice formation
     * A weighted average of the original damage, weighted by the old concentration, and the
     * undamaged new ice damage (1), weighted by the concentration of new ice. */
    damage[i] = (damage[i] * (cice[i] - lateralGrowth) + lateralGrowth) / cice[i];

    /* 2. Constant healing
     * Damage healing using a constant timescale. Originally conceived as an exponential decay, but
     * then revised to a linear one. */
    // This is what Sylvain and Pierre (Bouillon and Rampal, 2015)
    // damage[i] +=  damage[i] * tstep.step / tD;

    // This is what Véro did (Dansereau et al., 2016)
    damage[i] +=  tstep.step / tD;
    damage[i] = std::min(1., damage[i]);
}

}
