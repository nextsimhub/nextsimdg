/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConstantAtmosphereBoundary.hpp"

namespace Nextsim {

ConstantAtmosphereBoundary::ConstantAtmosphereBoundary()
    : IAtmosphereBoundary()
{
}

void ConstantAtmosphereBoundary::setData(const ModelState::DataMap& ms)
{
    // Directly set the array values
    qiaAccessor.getHostRW() = 305.288; // Pulled from NSColumnPhysics_test.cpp: New Ice Formation
    dqia_dtAccessor.getHostRW() = 4.5036;
    qowAccessor.getHostRW() = 307.546;
    sublAccessor.getHostRW() = 0.; // Seems unlikely…
    snowAccessor.getHostRW() = 0.;
    rainAccessor.getHostRW() = 0.;
    evapAccessor.getHostRW() = 0; // somehow...
    uwindAccessor.getHostRW() = 0;
    vwindAccessor.getHostRW() = 0;
    penSWAccessor.getHostRW() = 0.;
}

void ConstantAtmosphereBoundary::update(const TimestepTime& tst)
{
    IAtmosphereBoundary::update(tst);
}

} /* namespace Nextsim */
