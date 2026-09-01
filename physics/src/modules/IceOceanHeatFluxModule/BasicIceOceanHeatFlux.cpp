/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/BasicIceOceanHeatFlux.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

KERNEL_IMPL_FUNCTION static FloatType doOne(
    FloatType tBot, FloatType sst, FloatType mlBulkCp, FloatType timeT)
{
    // The ice-ocean flux is extremely sensitive to the small difference (sst - tf) and the large
    // multiplication by mlBulkCp / timeT amplifies any rounding error. In single precision this
    // difference suffers from catastrophic cancellation, which can drive spurious freezing/melting
    // and cause runaway ice growth in longer simulations.
    const double sstD = static_cast<double>(sst);
    const double tBotD = static_cast<double>(tBot);
    const double mlBulkCpD = static_cast<double>(mlBulkCp);
    const double timeTD = static_cast<double>(timeT);
    // Transfer rate depends on the mixed layer depth and the relaxation time scale
    return static_cast<FloatType>((sstD - tBotD) * mlBulkCpD / timeTD);
}

void BasicIceOceanHeatFlux::update(const TimestepTime& tst)
{
    auto& qio = qioAccessor.getAutoRW();
    const auto& cice = ciceAccessor.getAutoRO();
    const auto& mlBulkCp = mlBulkCpAccessor.getAutoRO();
    const auto& sst = sstAccessor.getAutoRO();
    const auto& tf = tfAccessor.getAutoRO();

    const FloatType dt = tst.step.seconds();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // Use the timestep length as the relaxation time scale
        if (cice[i] > 0.0_ft) {
            qio[i] = doOne(tf[i], sst[i], mlBulkCp[i], dt);
        } else {
            qio[i] = 0.0_ft;
        }
    });
}

} /* namespace Nextsim */
