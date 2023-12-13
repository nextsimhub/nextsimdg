/*!
 * @file ConstantOceanBoundary.cpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConstantOceanBoundary.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/Module.hpp"
#include "include/constants.hpp"

namespace Nextsim {
ConstantOceanBoundary::ConstantOceanBoundary()
    : IOceanBoundary()
{
}

void ConstantOceanBoundary::setData(const ModelState::DataMap& ms)
{
    IOceanBoundary::setData(ms);
    // Directly set the array values
    sss = 32.;
    u = 0;
    v = 0;
    mld = 10.;
    double tf32 = -1.751; // Hand calculated from S = 32 using UNESCO
    tf = tf32;
    sst = tf32; // Tf == SST ensures that there is no ice-ocean heat flux
    cpml = Water::cp * Water::rho * mld;
    qio = 0.;
}

void ConstantOceanBoundary::updateBefore(const TimestepTime& tst)
{
    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void ConstantOceanBoundary::updateAfter(const TimestepTime& tst) { }
} /* namespace Nextsim */
