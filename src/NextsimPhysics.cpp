/*!
 * @file NextsimPhysics.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NextsimPhysics.hpp"

#include <cmath>
#include "include/PrognosticData.hpp"
#include "include/ExternalData.hpp"
#include "include/PhysicsData.hpp"
#include "include/ElementData.hpp"

#include "include/constants.hpp"

namespace Nextsim {

NextsimPhysics::SpecificHumidity NextsimPhysics::specificHumidityWater;
NextsimPhysics::SpecificHumidityIce NextsimPhysics::specificHumidityIce;
double NextsimPhysics::dragOcean_q;

double stefanBoltzmannLaw(double temperature);

//NextsimPhysics::NextsimPhysics()
//{
//    // TODO Auto-generated constructor stub
//
//}
//
//NextsimPhysics::~NextsimPhysics()
//{
//    // TODO Auto-generated destructor stub
//}

void updateDerivedDataImpl(
        const PrognosticData& prog,
        const ExternalData& exter,
        PhysicsData& phys)
{
    phys.specificHumidityAir() = NextsimPhysics::specificHumidityWater(exter.airTemperature(), exter.airPressure());
    phys.specificHumidityWater() = NextsimPhysics::specificHumidityWater(prog.seaSurfaceTemperature(), exter.airPressure(), prog.seaSurfaceSalinity());

    phys.airDensity() = exter.airPressure() / (Air::Ra * kelvin(exter.airTemperature()));
}
void NextsimPhysics::updateDerivedDataStatic(ElementData& data)
{
    updateDerivedDataImpl(data, data, data);
}

void massFluxOpenWaterImpl(const PrognosticData& prog, PhysicsData& phys)
{
    double specificHumidityDifference = phys.specificHumidityWater() - phys.specificHumidityAir();
    phys.evaporationRate() = NextsimPhysics::dragOcean_q * phys.airDensity() * phys.windSpeed()
            * specificHumidityDifference;
}
void NextsimPhysics::massFluxOpenWaterStatic(ElementData& data)
{
    massFluxOpenWaterImpl(data, data);
}

void momentumFluxOpenWaterImpl(const PrognosticData& prog, PhysicsData& phys)
{
    phys.dragPressure() = phys.airDensity() * NextsimPhysics::dragOcean_m(phys.windSpeed());
}
void NextsimPhysics::momentumFluxOpenWaterStatic(ElementData& data)
{
    momentumFluxOpenWaterImpl(data, data);
}

void heatFluxOpenWaterImpl(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys)
{
    // Latent heat flux from evaporation and latent heat
    phys.QLatentHeat() = phys.evaporationRate() * NextsimPhysics::latentHeatWater(prog.seaSurfaceTemperature());

    // Sensible heat
    double heatCapacity = Air::cp + phys.specificHumidityAir()*Water::cp;
    phys.QSensibleHeat() = NextsimPhysics::dragOcean_t * phys.airDensity() * heatCapacity * phys.windSpeed() *
            (prog.seaSurfaceTemperature() - exter.airTemperature());

    // Shortwave flux
    phys.QShortwave() = -exter.incomingShortwave() * (1 - phys.oceanAlbedo());

    // Longwave flux
    phys.QLongwave() = stefanBoltzmannLaw(prog.seaSurfaceTemperature()) - exter.incomingLongwave();

    // Total flux
    phys.QOpenWater() = phys.QLatentHeat() +
            phys.QSensibleHeat() +
            phys.QLongwave() +
            phys.QShortwave();
}
void NextsimPhysics::heatFluxOpenWaterStatic(ElementData& data)
{
    heatFluxOpenWaterImpl(data, data, data);
}

void massFluxIceAtmosphere(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys)
{

}

void heatFluxIceAtmosphere(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys)
{

}

void massFluxIceOcean(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys)
{

}

void heatFluxIceOcean(const PrognosticData& prog, const ExternalData &exter, PhysicsData& phys)
{

}

void NextsimPhysics::setDragOcean_q(double doq)
{
    dragOcean_q = doq;
}

void NextsimPhysics::setDragOcean_t(double dot)
{
    dragOcean_t = dot;
}

double NextsimPhysics::dragOcean_m(double windSpeed)
{
    // Drag coefficient from Gill(1982) / Smith (1980)
    return 1e-3 * std::max(1.0, std::min(2.0, 0.61 + 0.063 * windSpeed));
}

double NextsimPhysics::latentHeatWater(double temperature)
{
    // Polynomial approximation expressed using Horner's scheme
    return Water::Lv0 +
            temperature * (-2.36418e3 +
                    temperature * (1.58927 +
                            temperature * (- 6.14342e-2)));
}

NextsimPhysics::SpecificHumidity::SpecificHumidity()
    : SpecificHumidity(
            6.1121e2, 18.729, 257.87, 227.3,
            7.2e-4, 3.20e-6, 5.9e-10)
{ }

// alpha (and beta) are consistent between both ice and water
NextsimPhysics::SpecificHumidity::SpecificHumidity(
        double a, double b, double c, double d,
        double bigA, double bigB, double bigC)
    : m_a(a)
    , m_b(b)
    , m_c(c)
    , m_d(d)
    , m_bigA(bigA)
    , m_bigB(bigB)
    , m_bigC(bigC)
    , m_alpha(0.62197)
    , m_beta(1 - m_alpha)
{ }

double NextsimPhysics::SpecificHumidity::operator()(const double temperature, const double pressure) const
{
    return this->operator ()(temperature, 0); // Zero salinity
}

double NextsimPhysics::SpecificHumidity::operator()(const double temperature, const double pressure, const double salinity) const
{
    double estCalc = est(temperature, salinity);
    double fCalc = f(temperature, pressure);
    double sphum = m_alpha * fCalc * estCalc / (pressure - m_beta * fCalc * estCalc);

    return sphum;
}

NextsimPhysics::SpecificHumidityIce::SpecificHumidityIce()
    : SpecificHumidity(
    6.1115e2, 23.036, 279.82, 333.7,
    2.2e-4, 3.83e-6, 6.4e-10)
{ }

double NextsimPhysics::SpecificHumidityIce::operator()(const double temperature, const double pressure) const
{
    return this->SpecificHumidity::operator()(temperature, pressure);
}

double NextsimPhysics::SpecificHumidityIce::dq_dT(const double temperature, const double pressure) const
{
    double df_dT = 2 * m_bigC * m_bigB * temperature;
    double numerator = m_a * m_b * m_c - temperature * (2 * m_c + temperature);
    double denominator = m_d * pow(m_c + temperature, 2);
    double estCalc = est(temperature, 0);
    double fCalc = f(pressure, temperature);
    double dest_dT = numerator / denominator * estCalc;
    numerator = m_alpha * pressure * (fCalc * dest_dT + estCalc * df_dT);
    denominator = pow(pressure - m_beta * estCalc * fCalc, 2);
    return numerator / denominator;
}

// Specific humidity terms
double NextsimPhysics::SpecificHumidity::f(const double pressurePa, const double temperature) const
{
    double pressure_mb = pressurePa * 0.01;
    return 1 + m_bigA + pressure_mb * (m_bigB + m_bigC * temperature * temperature);
}

double NextsimPhysics::SpecificHumidity::est(const double temperature, const double salinity) const
{
    double salFactor = 1 - 5.37e-4 * salinity;
    return m_a * exp( (m_b - temperature / m_d) * temperature / (temperature + m_c) ) * salFactor;
}

double stefanBoltzmannLaw(double temperatureC)
{
    return PhysicalConstants::sigma * std::pow(kelvin(temperatureC), 4);
}
} /* namespace Nextsim */
