/*!
 * @file BasicIceOceanHeatFlux.cpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BasicIceOceanHeatFlux.hpp"

#include "include/ExternalData.hpp"
#include "include/PrognosticData.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double BasicIceOceanHeatFlux::flux(
    const PrognosticData& prog, const ExternalData& exter, const PhysicsData& phys, const NextsimPhysics& nsp)
{
    // The ice bottom temperature is the freezing point of the surface seawater
    double iceBottomTemperature = prog.freezingPoint();
    double tDiff = prog.seaSurfaceTemperature() - iceBottomTemperature;

    // Transfer rate depends on the mixed layer depth an d a timescale. Here, it is the timestep
    return tDiff * exter.mixedLayerBulkHeatCapacity() / prog.timestep();
}

} /* namespace Nextsim */
