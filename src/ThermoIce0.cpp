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
    // True constants
    const double freezingPointIce = -Water::mu * Ice::s;
    const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
    const double bulkLHFusionIce = Water::Lf * Ice::rho;

    // TODO: Move these to the runtime constant class
    const double k_s = 0.3096; // Thermal conductivity of snow
    const bool doFlooding = true; // Whether flooded snow converts to ice

    if (prog.iceThickness() == 0 || prog.iceConcentration() == 0) {
        phys.iceTrueThickness() = 0;
        phys.snowTrueThickness() = 0;
        phys.updatedIceSurfaceTemperature() = freezingPointIce;

        return;
    }

    // Calculate the true slab thickness from the effective thickness
    phys.iceTrueThickness() = prog.iceThickness() / prog.iceConcentration();
    phys.snowTrueThickness() = prog.snowThickness() / prog.iceConcentration();

    double oldIceThickness = phys.iceTrueThickness();

    double iceTemperature = prog.iceTemperatures()[0];
    // Heat transfer coefficient
    double k_lSlab = k_s * Ice::kappa /
            (k_s * phys.iceTrueThickness() + Ice::kappa * phys.snowTrueThickness());
    double QIceConduction = k_lSlab * (exter.iceBottomTemperature() - iceTemperature);
    double remainingFlux = QIceConduction - phys.QIceAtmosphere();
    phys.updatedIceSurfaceTemperature() = iceTemperature +
            remainingFlux /
            (k_lSlab + phys.QDerivativeWRTTemperature());

    // Clamp the maximum temperature of the ice to the melting point of ice or snow
    double meltingLimit = (phys.snowTrueThickness() > 0.) ? 0 : freezingPointIce;
    phys.updatedIceSurfaceTemperature() =
            std::min(meltingLimit, phys.updatedIceSurfaceTemperature());

    // Top melt. Melting rate is non-positive.
    double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow; // [m³ s⁻¹]
    double snowSublRate = phys.sublimationRate() / Ice::rhoSnow; // [m³ s⁻¹]

    phys.snowTrueThickness() += (snowMeltRate - snowSublRate) * prog.timestep();
    // Use excess flux to melt ice. Non-positive value
    double excessIceMelt = std::min(phys.snowTrueThickness(), 0.) *
            bulkLHFusionSnow/bulkLHFusionIce;
    // With the excess flux noted, clamp the snow thickness to a minimum of zero.
    phys.snowTrueThickness() = std::max(phys.snowTrueThickness(), 0.);
    // Then add snowfall back on top
    phys.snowTrueThickness() += exter.snowfall() * prog.timestep() / Ice::rhoSnow;

    // Bottom melt or growth
    double iceBottomChange = (QIceConduction - phys.QIceOceanHeat()) * prog.timestep() /
            bulkLHFusionIce;

    // Total thickness change
    double iceThicknessChange = excessIceMelt + iceBottomChange;
    phys.iceTrueThickness() += iceThicknessChange;

    // Amount of melting (only) at the top and bottom of the ice
    double topMelt = std::min(excessIceMelt, 0.);
    double botMelt = std::min(iceBottomChange, 0.);

    // Snow to ice conversion
    double iceDraught = (phys.iceTrueThickness() * Ice::rho + phys.snowTrueThickness() * Ice::rhoSnow) /
            Water::rhoOcean;
    if (doFlooding && iceDraught > phys.iceTrueThickness()) {
        // Keep a running total of the ice formed from flooded snow
        double newIce = iceDraught - phys.iceTrueThickness();
        phys.totalIceFromSnow() += newIce;

        // Convert all the submerged snow to ice
        phys.iceTrueThickness() = iceDraught;
        phys.snowTrueThickness() -= newIce * Ice::rho / Ice::rhoSnow;
    }

    if (phys.iceTrueThickness() < ModelConstants::hMin) {
        // Reduce the melting to reach zero thickness, while keeping the
        // between top and bottom melting
        if (iceThicknessChange < 0) {
            double scaling = -oldIceThickness / iceThicknessChange;
            topMelt *= scaling;
            botMelt *= scaling;
        }

        // No snow was converted to ice
        phys.totalIceFromSnow() = 0;

        // Change in thickness is all of the old thickness
        iceThicknessChange = -oldIceThickness;

        // The ice-ocean flux includes all the latent heat
        phys.QIceOceanHeat() += phys.iceTrueThickness() * bulkLHFusionIce / prog.timestep()
                + phys.snowTrueThickness() * bulkLHFusionSnow / prog.timestep();

        // No ice, no snow and the surface temperature is the melting point of ice
        phys.iceTrueThickness() = 0;
        phys.snowTrueThickness() = 0;
        phys.updatedIceSurfaceTemperature() = freezingPointIce;
    }
}

} /* namespace Nextsim */
