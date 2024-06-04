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

template <>
const std::map<int, std::string> Configured<ConstantHealing>::keyMap
    = { { ConstantHealing::TD_KEY, "ConstantHealing.td" } };

void ConstantHealing::configure()
{
    tD = Configured::getConfiguration(keyMap.at(TD_KEY), tDDefault);
}

ModelState ConstantHealing::getStateRecursive(const Nextsim::OutputSpec& os) const
{
    return { {}, { { keyMap.at(TD_KEY), tD } } };
}

ConstantHealing::HelpMap& ConstantHealing::getHelpText(HelpMap& map, bool getAll)
{
    map["ConstantHealing"]
        = { { keyMap.at(TD_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(tDDefault),
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
void ConstantHealing::update(const TimestepTime& tstep) { }

}