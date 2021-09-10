/*!
 * @file NextsimPhysics.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NextsimPhysics.hpp"

#include <cmath>
#include "include/PrognosticData.hpp"
#include "include/PhysicsData.hpp"

namespace Nextsim {

NextsimPhysics::NextsimPhysics()
    : specificHumidityWater()
    , specificHumidityIce()
{
    // TODO Auto-generated constructor stub

}

NextsimPhysics::~NextsimPhysics()
{
    // TODO Auto-generated destructor stub
}

void NextsimPhysics::updateDerivedData(const PrognosticData& prog, PhysicsData& phys)
{

}

void NextsimPhysics::massFluxOpenWater(const PrognosticData& prog, PhysicsData& phys)
{
    double specificHumidityDifference = phys.specificHumidityWater() - phys.specificHumidityAir();
    phys.evaporationRate() = phys.dragOcean_q() * phys.airDensity() * phys.windSpeed()
            * specificHumidityDifference;
}

void NextsimPhysics::momentumFluxOpenWater(const PrognosticData& prog, PhysicsData& phys)
{
    phys.dragPressure() = phys.airDensity() * phys.dragOcean_m();
}

void NextsimPhysics::heatFluxOpenWater(const PrognosticData& prog, PhysicsData& phys)
{

}

NextsimPhysics::SpecificHumidity::SpecificHumidity()
    : m_a(6.1121e2)
    , m_b(18.729)
    , m_c(257.87)
    , m_d(227.3)
    , m_bigA(7.2e-4)
    , m_bigB(3.20e-6)
    , m_bigC(5.9e-10)
    , alpha(0.62197)
    , beta(1-alpha)
{ }

double NextsimPhysics::SpecificHumidity::operator()(double temperature, double pressure)
{
    return this->operator ()(temperature, 0); // Zero salinity
}

double NextsimPhysics::SpecificHumidity::operator()(double temperature, double pressure, double salinity)
{
    double p_mb = pressure * 0.01; // Convert from kPa to mbar
    double f = 1 + m_bigA + p_mb * (m_bigB + m_bigC * temperature * temperature);
    double est = m_a * exp( (m_b - temperature / m_d) * temperature / (temperature + m_c) ) *
            ( 1 - 5.37e-4 * salinity);
    double sphum = alpha * f * est / (pressure - beta * f * est);

    return sphum;
}

NextsimPhysics::SpecificHumidityIce::SpecificHumidityIce()
    : SpecificHumidity()
    , m_a(6.1115e2)
    , m_b(23.036)
    , m_c(279.82)
    , m_d(333.7)
    , m_bigA(2.2e-4)
    , m_bigB(3.83e-6)
    , m_bigC(6.4e-10)
{ }

double NextsimPhysics::SpecificHumidityIce::operator()(double temperature, double pressure)
{
    return this->SpecificHumidity::operator()(temperature, pressure);
}
} /* namespace Nextsim */
