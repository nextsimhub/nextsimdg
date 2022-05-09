/*!
 * @file AtmosphereOceanState.cpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/AtmosphereOceanState.hpp"

namespace Nextsim {

AtmosphereOceanState::AtmosphereOceanState()
{
    registerProtectedArray(ProtectedArray::HTRUE_ICE, &hTrueIce);
    registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hTrueSnow);
}

void AtmosphereOceanState::setData(const ModelState&) { }
ModelState AtmosphereOceanState::getState() const
{
    return {
        { "True ice thickness", hTrueIce },
        { "True snow thickness", hTrueSnow },
    };
}
ModelState AtmosphereOceanState::getState(const OutputLevel&) const { return getState(); }
std::string AtmosphereOceanState::getName() const { return "AtmosphereOceanState"; }
std::set<std::string> AtmosphereOceanState::hFields() const
{
    return { "True ice thickness", "True snow thickness" };
}

void AtmosphereOceanState::update(const TimestepTime& tst)
{
    hTrueSnow = hSnowCell / cIce;
    hTrueIce = hIceCell / cIce;
}

} /* namespace Nextsim */
