/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/BasicIceOceanHeatFlux.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

KERNEL_IMPL_FUNCTION static double doOne(double tBot, double sst, double mlBulkCp, double timeT)
{
    // Transfer rate depends on the mixed layer depth and the relaxation time scale
    return (sst - tBot) * mlBulkCp / timeT;
}

void BasicIceOceanHeatFlux::update(const TimestepTime& tst)
{
    auto& qio = qioAccessor.getAutoRW();
    const auto& cice = ciceAccessor.getAutoRO();
    const auto& mlBulkCp = mlBulkCpAccessor.getAutoRO();
    const auto& sst = sstAccessor.getAutoRO();
    const auto& tf = tfAccessor.getAutoRO();

    const double dt = tst.step.seconds();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Use the timestep length as the relaxation time scale
        if (cice[i] > 0.) {
            qio[i] = doOne(tf[i], sst[i], mlBulkCp[i], dt);
        } else {
            qio[i] = 0.;
        }
    });
}

} /* namespace Nextsim */
