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
    registerProtectedArray(ProtectedArray::HTRUE_ICE, &hTrueIce);
    registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hTrueSnow);
    atmosStateImpl->setData(ms);
    oceanStateImpl->setData(ms);
}
ModelState AtmosphereOceanState::getState() const
{
    return { {
                 { "True ice thickness", hTrueIce },
                 { "True snow thickness", hTrueSnow },
             },
        {} };
}
ModelState AtmosphereOceanState::getState(const OutputLevel&) const { return getState(); }
ModelState AtmosphereOceanState::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    state.merge(atmosStateImpl->getStateRecursive(os));
    state.merge(oceanStateImpl->getStateRecursive(os));
    return state;
}

std::string AtmosphereOceanState::getName() const { return "AtmosphereOceanState"; }
std::unordered_set<std::string> AtmosphereOceanState::hFields() const
{
    return { "True ice thickness", "True snow thickness" };
}

void AtmosphereOceanState::configure()
{
    atmosStateImpl = std::move(Module::getInstance<AtmosphereState>());
    tryConfigure(*atmosStateImpl);

    oceanStateImpl = std::move(Module::getInstance<OceanState>());
    tryConfigure(*oceanStateImpl);
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
    oceanStateImpl->update(tst);

    hTrueSnow = hSnowCell / cIce;
    hTrueIce = hIceCell / cIce;
}

} /* namespace Nextsim */
