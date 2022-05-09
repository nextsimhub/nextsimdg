/*!
 * @file BasicIceOceanHeatFlux.cpp
 *
 * @date Oct 19, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BasicIceOceanHeatFluxED.hpp"
#include "include/ExternalData.hpp"
#include "include/PrognosticElementData.hpp"
#include "include/constants.hpp"

namespace Nextsim {

static inline double doOne(double tBot, double sst, double mlBulkCp, double ts)
{
    // Transfer rate depends on the mixed layer depth an d a timescale. Here, it is the timestep
    return (sst - tBot) * mlBulkCp / ts;
}

double BasicIceOceanHeatFluxED::flux(const PrognosticElementData& prog, const ExternalData& exter,
    const PhysicsData& phys, const NextsimPhysics& nsp)
{
    // The ice bottom temperature is the freezing point of the surface seawater
    return doOne(prog.freezingPoint(), prog.seaSurfaceTemperature(),
        exter.mixedLayerBulkHeatCapacity(), prog.timestep());
}

} /* namespace Nextsim */
