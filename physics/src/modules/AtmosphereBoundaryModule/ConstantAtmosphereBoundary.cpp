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
    qia = 305.288; // Pulled from NSColumnPhysics_test.cpp: New Ice Formation
    dqia_dt = 4.5036;
    qow = 307.546;
    subl = 0.; // Seems unlikely…
    snow = 0.;
    rain = 0.;
    evap = 0; // somehow...
    uwind = 0;
    vwind = 0;
    penSW = 0.;
}

void ConstantAtmosphereBoundary::update(const TimestepTime& tst)
{
    IAtmosphereBoundary::update(tst);
}

} /* namespace Nextsim */
