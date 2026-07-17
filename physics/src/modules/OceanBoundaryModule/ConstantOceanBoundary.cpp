/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConstantOceanBoundary.hpp"

#include "include/Finalizer.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {
ConstantOceanBoundary::ConstantOceanBoundary()
    : IOceanBoundary()
{
}

void ConstantOceanBoundary::setData(const ModelState::DataMap& ms)
{
    Finalizer::registerUnique(Module::finalize<IIceOceanHeatFlux>);

    IOceanBoundary::setData(ms);
    // Directly set the array values
    sssAccessor.getHostRW() = 32.;
    uAccessor.getHostRW() = 0;
    vAccessor.getHostRW() = 0;
    HField& mld = mldAccessor.getHostRW();
    mld = 10.;
    FloatType tf32 = -1.751; // Hand calculated from S = 32 using UNESCO
    tfAccessor.getHostRW() = tf32;
    sstAccessor.getHostRW() = tf32; // Tf == SST ensures that there is no ice-ocean heat flux
    cpmlAccessor.getHostRW() = Water::cp * Water::rho * mld;
    qioAccessor.getHostRW() = 0.;
    sshAccessor.getHostRW() = 0.;
}

void ConstantOceanBoundary::updateBefore(const TimestepTime& tst)
{
    Module::getImplementation<IIceOceanHeatFlux>().update(tst);
}

void ConstantOceanBoundary::updateAfter(const TimestepTime& tst) { }
} /* namespace Nextsim */
