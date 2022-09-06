/*!
 * @file AtmosphereOceanState.cpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/AtmosphereOceanState.hpp"
#include "include/AtmosphereStateModule.hpp"
#include "include/OceanStateModule.hpp"

namespace Nextsim {

AtmosphereOceanState::AtmosphereOceanState() { }

void AtmosphereOceanState::setData(const ModelState::DataMap& ms)
{
    atmosStateImpl->setData(ms);
}
ModelState AtmosphereOceanState::getState() const { return { {}, {} }; }
ModelState AtmosphereOceanState::getState(const OutputLevel&) const { return getState(); }
ModelState AtmosphereOceanState::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    state.merge(atmosStateImpl->getStateRecursive(os));
    return state;
}

std::string AtmosphereOceanState::getName() const { return "AtmosphereOceanState"; }
std::unordered_set<std::string> AtmosphereOceanState::hFields() const
{
    return { };
}

void AtmosphereOceanState::configure()
{
    atmosStateImpl = std::move(Module::getInstance<AtmosphereState>());
    tryConfigure(*atmosStateImpl);
}

AtmosphereOceanState::HelpMap& AtmosphereOceanState::getHelpRecursive(HelpMap& map, bool getAll)
{
    Module::getHelpRecursive<AtmosphereState>(map, getAll);
    Module::getHelpRecursive<OceanState>(map, getAll);
    return map;
}

void AtmosphereOceanState::update(const TimestepTime& tst)
{
    atmosStateImpl->update(tst);
}

} /* namespace Nextsim */
