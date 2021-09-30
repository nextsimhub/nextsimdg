/*
 * @file ThermoIce0.cpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0.hpp"

#include "include/PrognosticData.hpp"
#include "include/ExternalData.hpp"
#include "include/PhysicsData.hpp"
#include "include/NextsimPhysics.hpp"

#include "include/constants.hpp"

namespace Nextsim {

void ThermoIce0::calculate(
        const PrognosticData& prog,
        const ExternalData& exter,
        PhysicsData& phys,
        NextsimPhysics& nsphys)
{
    const double freezingPointIce = -Water::mu * Ice::s;
    // TODO: Move this to the runtime constant class
    const double k_s = 0.3096; // Thermal conductivity of snow

    if (prog.iceThickness() == 0 || prog.iceConcentration() == 0) {
        phys.iceTrueThickness() = 0;
        phys.snowTrueThickness() = 0;
        phys.updatedIceSurfaceTemperature() = freezingPointIce;
    } else {
        // Calculate the true slab thickness from the effective thickness
        phys.iceTrueThickness() = prog.iceThickness() / prog.iceConcentration();
        phys.snowTrueThickness() = prog.snowThickness() / prog.iceConcentration();
    }
    // Heat transfer coefficient
    double k_lSlab = k_s * Ice::kappa /
            (k_s * phys.iceTrueThickness() + Ice::kappa * phys.snowTrueThickness());

}

} /* namespace Nextsim */
