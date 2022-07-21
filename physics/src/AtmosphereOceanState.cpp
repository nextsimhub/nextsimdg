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

void AtmosphereOceanState::setData(const ModelState& ms)
{
    atmosStateImpl->setData(ms);
    oceanStateImpl->setData(ms);
}
ModelState AtmosphereOceanState::getState() const { return ModelState(); }
ModelState AtmosphereOceanState::getState(const OutputLevel&) const { return getState(); }
std::string AtmosphereOceanState::getName() const { return "AtmosphereOceanState"; }
std::unordered_set<std::string> AtmosphereOceanState::hFields() const
{
    return std::unordered_set<std::string>();
}

void AtmosphereOceanState::configure()
{
    atmosStateImpl = std::move(Module::getInstance<AtmosphereState>());
    tryConfigure(*atmosStateImpl);

    oceanStateImpl = std::move(Module::getInstance<OceanState>());
    tryConfigure(*oceanStateImpl);
}

void AtmosphereOceanState::update(const TimestepTime& tst)
{
    atmosStateImpl->updateAtmos(tst);
    oceanStateImpl->updateOcean(tst);
}

} /* namespace Nextsim */
