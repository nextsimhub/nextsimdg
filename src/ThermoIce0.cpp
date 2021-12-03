/*!
 * @file ThermoIce0.cpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0.hpp"

#include "include/ExternalData.hpp"
#include "include/NextsimPhysics.hpp"
#include "include/PhysicsData.hpp"
#include "include/PrognosticData.hpp"

#include "include/constants.hpp"

namespace Nextsim {

double ThermoIce0::k_s = 0;
bool ThermoIce0::doFlooding = true;

template <>
const std::map<int, std::string> Configured<ThermoIce0>::keyMap = {
    { ThermoIce0::KS_KEY, "thermoice0.ks" },
    { ThermoIce0::FLOODING_KEY, "thermoice0.flooding" },
};

void ThermoIce0::configure()
{
    k_s = Configured::getConfiguration(keyMap.at(KS_KEY), 0.3096);
    doFlooding = Configured::getConfiguration(keyMap.at(FLOODING_KEY), true);
}

void ThermoIce0::calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys,
    NextsimPhysics& nsphys)
{
    // True constants
    const double freezingPointIce = -Water::mu * Ice::s;
    const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
    const double bulkLHFusionIce = Water::Lf * Ice::rho;

    // Initialize the updated snow thickness
    // phys.updatedSnowTrueThickness() = prog.snowTrueThickness();

    if (prog.iceThickness() == 0 || prog.iceConcentration() == 0) {
        phys.updatedIceTrueThickness() = 0;
        phys.updatedSnowTrueThickness() = 0;
        phys.updatedIceSurfaceTemperature() = freezingPointIce;

        return;
    }

    double oldIceThickness = prog.iceTrueThickness();

    double iceTemperature = prog.iceTemperatures()[0];
    double tBot = prog.freezingPoint();
    // Heat transfer coefficient
    double k_lSlab = k_s * Ice::kappa
        / (k_s * prog.iceTrueThickness() + Ice::kappa * prog.snowTrueThickness());
    double QIceConduction = k_lSlab * (tBot - iceTemperature);
    double remainingFlux = QIceConduction - nsphys.QIceAtmosphere();
    phys.updatedIceSurfaceTemperature()
        = iceTemperature + remainingFlux / (k_lSlab + nsphys.QDerivativeWRTTemperature());

    // Clamp the maximum temperature of the ice to the melting point of ice or snow
    double meltingLimit = (prog.snowTrueThickness() > 0.) ? 0 : freezingPointIce;
    phys.updatedIceSurfaceTemperature()
        = std::min(meltingLimit, phys.updatedIceSurfaceTemperature());

    // Top melt. Melting rate is non-positive.
    double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow; // [m³ s⁻¹]
    double snowSublRate = nsphys.sublimationRate() / Ice::rhoSnow; // [m³ s⁻¹]

    phys.updatedSnowTrueThickness() += (snowMeltRate - snowSublRate) * prog.timestep();
    // Use excess flux to melt ice. Non-positive value
    double excessIceMelt
        = std::min(phys.updatedSnowTrueThickness(), 0.) * bulkLHFusionSnow / bulkLHFusionIce;
    // With the excess flux noted, clamp the snow thickness to a minimum of zero.
    phys.updatedSnowTrueThickness() = std::max(phys.updatedSnowTrueThickness(), 0.);
    // Then add snowfall back on top
    phys.updatedSnowTrueThickness() += exter.snowfall() * prog.timestep() / Ice::rhoSnow;

    // Bottom melt or growth
    double iceBottomChange
        = (QIceConduction - nsphys.QIceOceanHeat()) * prog.timestep() / bulkLHFusionIce;
    // Total thickness change
    double iceThicknessChange = excessIceMelt + iceBottomChange;
    phys.updatedIceTrueThickness() += iceThicknessChange;

    // Amount of melting (only) at the top and bottom of the ice
    double topMelt = std::min(excessIceMelt, 0.);
    double botMelt = std::min(iceBottomChange, 0.);

    // Snow to ice conversion
    double iceDraught = (phys.updatedIceTrueThickness() * Ice::rho
                            + phys.updatedSnowTrueThickness() * Ice::rhoSnow)
        / Water::rhoOcean;
    if (doFlooding && iceDraught > phys.updatedIceTrueThickness()) {
        // Keep a running total of the ice formed from flooded snow
        double newIce = iceDraught - phys.updatedIceTrueThickness();
        nsphys.totalIceFromSnow() += newIce;

        // Convert all the submerged snow to ice
        phys.updatedIceTrueThickness() = iceDraught;
        phys.updatedSnowTrueThickness() -= newIce * Ice::rho / Ice::rhoSnow;
    }

    if (phys.updatedIceTrueThickness() < NextsimPhysics::minimumIceThickness()) {
        // Reduce the melting to reach zero thickness, while keeping the
        // between top and bottom melting
        if (iceThicknessChange < 0) {
            double scaling = -oldIceThickness / iceThicknessChange;
            topMelt *= scaling;
            botMelt *= scaling;
        }

        // No snow was converted to ice
        nsphys.totalIceFromSnow() = 0;

        // Change in thickness is all of the old thickness
        iceThicknessChange = -oldIceThickness;

        // The ice-ocean flux includes all the latent heat
        nsphys.QIceOceanHeat() += phys.updatedIceTrueThickness() * bulkLHFusionIce / prog.timestep()
            + phys.updatedSnowTrueThickness() * bulkLHFusionSnow / prog.timestep();

        // No ice, no snow and the surface temperature is the melting point of ice
        phys.updatedIceTrueThickness() = 0;
        phys.updatedSnowTrueThickness() = 0;
        phys.updatedIceSurfaceTemperature() = freezingPointIce;
    }
}

} /* namespace Nextsim */
