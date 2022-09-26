/*!
 * @file ConstantAtmosphereBoundary.cpp
 *
 * @date Sep 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConstantAtmosphereBoundary.hpp"

namespace Nextsim {

ConstantAtmosphereBoundary::ConstantAtmosphereBoundary()
    : IAtmosphereBoundary()
{
}

void ConstantAtmosphereBoundary::setData(const ModelState::DataMap& ms)
{ /* nope! */
}

void ConstantAtmosphereBoundary::update(const TimestepTime& tst)
{
    IAtmosphereBoundary::update(tst);
    // Directly set the array values
    qia = 305.288; // Pulled from IceGrowth_test.cpp: New Ice Formation
    dqia_dt = 4.5036;
    sw_in = 0; // night
    lw_in = 311; // -1˚C air temperature
    precip = 0;
    evap = 0; // somehow...
    uwind = 0;
    vwind = 0;
}

} /* namespace Nextsim */
