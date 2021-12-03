/*!
 * @file NextsimPhysics.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NextsimPhysics.hpp"

#include "include/ElementData.hpp"
#include "include/ExternalData.hpp"
#include "include/PhysicsData.hpp"
#include "include/PrognosticData.hpp"
#include <cmath>

#include "include/IConcentrationModel.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/IThermodynamics.hpp"

#include "include/ModuleLoader.hpp"

#include "include/constants.hpp"

namespace Nextsim {

NextsimPhysics::SpecificHumidity NextsimPhysics::specHumWater;
NextsimPhysics::SpecificHumidityIce NextsimPhysics::specHumIce;
double NextsimPhysics::dragOcean_q;
double NextsimPhysics::dragOcean_t;
double NextsimPhysics::dragIce_t;
double NextsimPhysics::m_oceanAlbedo;
double NextsimPhysics::m_I0;
double NextsimPhysics::minc;
double NextsimPhysics::minh;

IIceOceanHeatFlux* NextsimPhysics::iceOceanHeatFluxImpl = nullptr;
IIceAlbedo* NextsimPhysics::iIceAlbedoImpl = nullptr;
IThermodynamics* NextsimPhysics::iThermo = nullptr;
IConcentrationModel* NextsimPhysics::iConcentrationModelImpl = nullptr;

double stefanBoltzmannLaw(double temperature);

NextsimPhysics::NextsimPhysics()
    : m_Qio(0)
    , m_newice(0)
{
}

template <>
const std::map<int, std::string> Configured<NextsimPhysics>::keyMap = {
    { NextsimPhysics::DRAGOCEANQ_KEY, "nextsim_thermo.drag_ocean_q" },
    { NextsimPhysics::DRAGOCEANT_KEY, "nextsim_thermo.drag_ocean_t" },
    { NextsimPhysics::DRAGICET_KEY, "nextsim_thermo.drag_ice_t" },
    { NextsimPhysics::OCEANALBEDO_KEY, "nextsim_thermo.albedoW" },
    { NextsimPhysics::I0_KEY, "nextsim_thermo.I_0" },
    { NextsimPhysics::MINC_KEY, "nextsim_thermo.min_conc" },
    { NextsimPhysics::MINH_KEY, "nextsim_thermo.min_thick" },
};

void NextsimPhysics::configure()
{
    ModuleLoader& loader = ModuleLoader::getLoader();

    iceOceanHeatFluxImpl = &loader.getImplementation<IIceOceanHeatFlux>();
    tryConfigure(iceOceanHeatFluxImpl);

    iIceAlbedoImpl = &loader.getImplementation<IIceAlbedo>();
    tryConfigure(iIceAlbedoImpl);

    iThermo = &loader.getImplementation<IThermodynamics>();
    tryConfigure(iThermo);

    iConcentrationModelImpl = &loader.getImplementation<IConcentrationModel>();
    tryConfigure(iConcentrationModelImpl);

    dragOcean_q = Configured::getConfiguration(keyMap.at(DRAGOCEANQ_KEY), 1.5e-3);
    dragOcean_t = Configured::getConfiguration(keyMap.at(DRAGOCEANT_KEY), 0.83e-3);
    dragIce_t = Configured::getConfiguration(keyMap.at(DRAGICET_KEY), 1.3e-3);
    m_oceanAlbedo = Configured::getConfiguration(keyMap.at(OCEANALBEDO_KEY), 0.07);
    m_I0 = Configured::getConfiguration(keyMap.at(I0_KEY), 0.17);
    minc = Configured::getConfiguration(keyMap.at(MINC_KEY), 1e-12);
    minh = Configured::getConfiguration(keyMap.at(MINH_KEY), 0.01);
}

void NextsimPhysics::updateSpecificHumidityAir(const ExternalData& exter, PhysicsData& phys)
{
    phys.specificHumidityAir() = specHumWater(exter.dewPoint2m(), exter.airPressure());
};

void NextsimPhysics::updateSpecificHumidityWater(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    phys.specificHumidityWater() = specHumWater(
        prog.seaSurfaceTemperature(), exter.airPressure(), prog.seaSurfaceSalinity());
}

void NextsimPhysics::updateSpecificHumidityIce(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    phys.specificHumidityIce() = specHumIce(prog.iceTemperatures()[0], exter.airPressure());
}

void NextsimPhysics::updateAirDensity(const ExternalData& exter, PhysicsData& phys)
{
    double Ra_wet = Air::Ra / (1 - phys.specificHumidityAir() * ( 1 - Vapour::Ra / Air::Ra));
    phys.airDensity() = exter.airPressure() / (Ra_wet * kelvin(exter.airTemperature()));
}

void NextsimPhysics::updateHeatCapacityWetAir(const ExternalData& exter, PhysicsData& phys)
{
    phys.heatCapacityWetAir() = Air::cp + phys.specificHumidityAir() * Water::cp;
};

void NextsimPhysics::calculate(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    massFluxOpenWater(phys);
    momentumFluxOpenWater(phys);
    heatFluxOpenWater(prog, exter, phys);

    massFluxIceAtmosphere(prog, phys);
    // Ice momentum fluxes are handled by the dynamics
    heatFluxIceAtmosphere(prog, exter, phys);

    // The mass flux is driven by the heat flux, so that is called first
    heatFluxIceOcean(prog, exter, phys);
    // Ice momentum fluxes are handled by the dynamics
    massFluxIceOcean(prog, exter, phys);
}

void NextsimPhysics::massFluxOpenWater(PhysicsData& phys)
{
    double specificHumidityDifference = phys.specificHumidityWater() - phys.specificHumidityAir();
    m_evap = dragOcean_q * phys.airDensity() * phys.windSpeed() * specificHumidityDifference;
}

void NextsimPhysics::momentumFluxOpenWater(PhysicsData& phys)
{
    phys.dragPressure() = phys.airDensity() * dragOcean_m(phys.windSpeed());
}

void NextsimPhysics::heatFluxOpenWater(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    // Latent heat flux from evaporation and condensation
    m_Qlhow = m_evap * latentHeatWater(prog.seaSurfaceTemperature());

    // Sensible heat flux
    m_Qshow = dragOcean_t * phys.airDensity() * phys.heatCapacityWetAir() * phys.windSpeed()
        * (prog.seaSurfaceTemperature() - exter.airTemperature());

    // Shortwave flux
    m_Qswow = -exter.incomingShortwave() * (1 - m_oceanAlbedo);

    // Longwave flux
    m_Qlwow = stefanBoltzmannLaw(prog.seaSurfaceTemperature()) - exter.incomingLongwave();

    // Total flux
    m_Qow = m_Qlhow + m_Qshow + m_Qlwow + m_Qswow;

}

void NextsimPhysics::massFluxIceAtmosphere(const PrognosticData& prog, PhysicsData& phys)
{
    m_subl = dragIce_t * phys.airDensity() * phys.windSpeed()
        * (phys.specificHumidityIce() - phys.specificHumidityAir());
}

void NextsimPhysics::heatFluxIceAtmosphere(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    // Latent heat flux from sublimation
    m_Qlhi = m_subl * latentHeatIce(prog.iceTemperatures()[0]);
    double dmdot_dT = dragIce_t * phys.airDensity() * phys.windSpeed()
        * specHumIce.dq_dT(prog.iceTemperatures()[0], exter.airPressure());
    m_dQ_dT += latentHeatIce(prog.iceTemperatures()[0]) * dmdot_dT;

    // Sensible heat flux
    m_Qshi = dragIce_t * phys.airDensity() * phys.heatCapacityWetAir() * phys.windSpeed()
        * (prog.iceTemperatures()[0] - exter.airTemperature());
    m_dQ_dT += dragIce_t * phys.airDensity() * phys.heatCapacityWetAir() * phys.windSpeed();

    // Shortwave flux
    double albedoValue = iIceAlbedoImpl->albedo(prog.iceTemperatures()[0],
        (prog.iceConcentration() > 0) ? (prog.snowThickness() / prog.iceConcentration()) : 0.);
    m_Qswi = -exter.incomingShortwave() * (1. - m_I0) * (1 - albedoValue);

    // Longwave flux
    m_Qlwi = stefanBoltzmannLaw(prog.iceTemperatures()[0]) - exter.incomingLongwave();
    m_dQ_dT += 4 / kelvin(prog.iceTemperatures()[0]) * stefanBoltzmannLaw(prog.iceTemperatures()[0]);

    // Total flux
    m_Qia = m_Qlhi + m_Qshi + m_Qlwi + m_Qswi;
}

void NextsimPhysics::massFluxIceOcean(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    m_hifroms = 0;

    iThermo->calculate(prog, exter, phys, *this);
    newIceFormation(prog, exter, phys);

    lateralGrowth(prog, exter, phys);

    // Apply the lower limit of concentration and thickness
    if (phys.updatedIceConcentration() < minc || phys.updatedIceTrueThickness() < minh) {
        m_Qow += phys.updatedIceConcentration() * Water::Lf
            * (phys.updatedIceTrueThickness() * Ice::rho
                + phys.updatedSnowTrueThickness() * Ice::rhoSnow)
            / prog.timestep();
        phys.updatedIceConcentration() = 0;
        phys.updatedIceTrueThickness() = 0;
        phys.updatedSnowTrueThickness() = 0;
    }
}

void NextsimPhysics::heatFluxIceOcean(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    m_Qio = iceOceanHeatFluxImpl->flux(prog, exter, phys, *this);
}

void NextsimPhysics::newIceFormation(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    // Flux cooling the ocean from open water
    // TODO Add assimilation fluxes here
    double coolingFlux = m_Qow;
    // Temperature change of the mixed layer during this timestep
    double deltaTml = -coolingFlux / exter.mixedLayerBulkHeatCapacity() * prog.timestep();
    // Initial temperature
    double t0 = prog.seaSurfaceTemperature();
    // Freezing point temperature
    double tf = prog.freezingPoint();
    // Final temperature
    double t1 = t0 + deltaTml;

    // deal with cooling below the freezing point
    if (t1 < tf) {
        // Heat lost cooling the mixed layer to freezing point
        double sensibleFlux = (tf - t0) / deltaTml * coolingFlux;
        // Any heat beyond that is latent heat forming new ice
        double latentFlux = coolingFlux - sensibleFlux;

        m_Qow = sensibleFlux;
        m_newice
            = latentFlux * prog.timestep() * (1 - prog.iceConcentration()) / (Ice::Lf * Ice::rho);
    }
}

// Update thickness with concentration
double updateThickness(double& thick, double oldConc, double deltaC, double deltaV)
{
    return thick += (deltaV - thick * deltaC) / (oldConc + deltaC);
}

void NextsimPhysics::lateralGrowth(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    NextsimPhysics& nsphys = *this;
    double del_c = 0; // Change in concentration due to lateral growth
    del_c += iConcentrationModelImpl->freeze(prog, phys, nsphys);
    if (phys.updatedIceTrueThickness() < prog.iceTrueThickness()) {
        del_c += iConcentrationModelImpl->melt(prog, phys, nsphys);
    }

    // Correct the ice thickness, snow thickness and open water flux based on the change in
    // concentration
    phys.updatedIceConcentration() = prog.iceConcentration() + del_c;

    if (phys.updatedIceConcentration() >= minc) {
        // The updated ice thickness must conserve volume
        updateThickness(phys.updatedIceTrueThickness(), prog.iceConcentration(), del_c, m_newice);

        if (del_c < 0) {
            // Snow is lost if the concentration decreases, and energy is returned to the ocean
            m_Qow -= del_c * phys.updatedSnowTrueThickness() * Water::Lf * Ice::rhoSnow
                / prog.timestep();
        } else {
            // Currently no new snow is implemented
            updateThickness(phys.updatedSnowTrueThickness(), prog.iceConcentration(), del_c, 0.);
        }
    }
}

double NextsimPhysics::dragOcean_m(double windSpeed)
{
    // Drag coefficient from Gill(1982) / Smith (1980)
    return 1e-3 * std::fmax(1.0, std::fmin(2.0, 0.61 + 0.063 * windSpeed));
}

double NextsimPhysics::latentHeatWater(double temperature)
{
    // Polynomial approximation expressed using Horner's scheme
    return Water::Lv0
        + temperature * (-2.36418e3 + temperature * (1.58927 + temperature * (-6.14342e-2)));
}

double NextsimPhysics::latentHeatIce(double temperature)
{
    return Water::Lv0 + Water::Lf - 240. + temperature * (-290. + temperature * (-4.));
}

NextsimPhysics::SpecificHumidity::SpecificHumidity()
    : SpecificHumidity(6.1121e2, 18.729, 257.87, 227.3, 7.2e-4, 3.20e-6, 5.9e-10)
{
}

// alpha (and beta) are consistent between both ice and water
NextsimPhysics::SpecificHumidity::SpecificHumidity(
    double a, double b, double c, double d, double bigA, double bigB, double bigC)
    : m_a(a)
    , m_b(b)
    , m_c(c)
    , m_d(d)
    , m_bigA(bigA)
    , m_bigB(bigB)
    , m_bigC(bigC)
    , m_alpha(0.62197)
    , m_beta(1 - m_alpha)
{
}

double NextsimPhysics::SpecificHumidity::operator()(
    const double temperature, const double pressure) const
{
    return this->operator()(temperature, pressure, 0); // Zero salinity
}

double NextsimPhysics::SpecificHumidity::operator()(
    const double temperature, const double pressure, const double salinity) const
{
    double estCalc = est(temperature, salinity);
    double fCalc = f(temperature, pressure);
    double sphum = m_alpha * fCalc * estCalc / (pressure - m_beta * fCalc * estCalc);

    return sphum;
}

NextsimPhysics::SpecificHumidityIce::SpecificHumidityIce()
    : SpecificHumidity(6.1115e2, 23.036, 279.82, 333.7, 2.2e-4, 3.83e-6, 6.4e-10)
{
}

double NextsimPhysics::SpecificHumidityIce::operator()(
    const double temperature, const double pressure) const
{
    return this->SpecificHumidity::operator()(temperature, pressure);
}

double NextsimPhysics::SpecificHumidityIce::dq_dT(
    const double temperature, const double pressure) const
{
    double df_dT = 2 * m_bigC * m_bigB * temperature;
    double numerator = m_b * m_c * m_d - temperature * (2 * m_c + temperature);
    double denominator = m_d * pow(m_c + temperature, 2);
    double estCalc = est(temperature, 0);
    double fCalc = f(temperature, pressure);
    double dest_dT = numerator / denominator * estCalc;
    numerator = m_alpha * pressure * (fCalc * dest_dT + estCalc * df_dT);
    denominator = pow(pressure - m_beta * estCalc * fCalc, 2);
    return numerator / denominator;
}

// Specific humidity terms
double NextsimPhysics::SpecificHumidity::f(const double temperature, const double pressurePa) const
{
    double pressure_mb = pressurePa * 0.01;
    return 1 + m_bigA + pressure_mb * (m_bigB + m_bigC * temperature * temperature);
}

double NextsimPhysics::SpecificHumidity::est(const double temperature, const double salinity) const
{
    double salFactor = 1 - 5.37e-4 * salinity;
    return m_a * exp((m_b - temperature / m_d) * temperature / (temperature + m_c)) * salFactor;
}

double stefanBoltzmannLaw(double temperatureC)
{
    return Ice::epsilon * PhysicalConstants::sigma * std::pow(kelvin(temperatureC), 4);
}
} /* namespace Nextsim */
