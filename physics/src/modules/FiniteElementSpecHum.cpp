/*
 * @file FiniteElementSpecHum.cpp
 *
 * @date May 3, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FiniteElementSpecHum.hpp"

#include <cmath>

namespace Nextsim {

FiniteElementSpecHum FiniteElementSpecHum::m_water(
    6.1121e2, 18.729, 257.87, 227.3, 7.2e-4, 3.20e-6, 5.9e-10);
FiniteElementSpecHum FiniteElementSpecHum::m_ice(
    6.1115e2, 23.036, 279.82, 333.7, 2.2e-4, 3.83e-6, 6.4e-10);

// Default construction constructs a water instance
FiniteElementSpecHum::FiniteElementSpecHum()
    : FiniteElementSpecHum(6.1121e2, 18.729, 257.87, 227.3, 7.2e-4, 3.20e-6, 5.9e-10)
{
}

// General constructor
FiniteElementSpecHum::FiniteElementSpecHum(
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

double FiniteElementSpecHum::operator()(double temperature, double pressure) const
{
    return operator()(temperature, pressure, 0);
}

double FiniteElementSpecHum::operator()(double temperature, double pressure, double salinity) const
{
    return calculate(temperature, pressure, salinity, false).first;
}

std::pair<double, double> FiniteElementSpecHum::valueAndDerivative(
    double temperature, double pressure) const
{
    return valueAndDerivative(temperature, pressure, 0);
}

std::pair<double, double> FiniteElementSpecHum::valueAndDerivative(
    double temperature, double pressure, double salinity) const
{
    return calculate(temperature, pressure, salinity, true);
}

std::pair<double, double> FiniteElementSpecHum::calculate(
    double temperature, double pressure, double salinity, bool doDeriv) const
{
    double estCalc = est(temperature, salinity);
    double fCalc = f(temperature, pressure);
    double sphum = m_alpha * fCalc * estCalc / (pressure - m_beta * fCalc * estCalc);

    double deriv = 0;

    if (doDeriv) {
        double df_dT = 2 * m_bigC * m_bigB * temperature;
        double numerator = m_b * m_c * m_d - temperature * (2 * m_c + temperature);
        double sqrtDenom = m_c + temperature;
        double denominator = m_d * sqrtDenom * sqrtDenom;
        double dest_dT = numerator / denominator * estCalc;
        numerator = m_alpha * pressure * (fCalc * dest_dT + estCalc * df_dT);
        sqrtDenom = pressure - m_beta * estCalc * fCalc;
        denominator = sqrtDenom * sqrtDenom;

        deriv = numerator / denominator;
    }

    return std::pair<double, double>(sphum, deriv);
}

// Specific humidity terms
double FiniteElementSpecHum::f(double temperature, double pressurePa) const
{
    double pressure_mb = pressurePa * 0.01;
    return 1 + m_bigA + pressure_mb * (m_bigB + m_bigC * temperature * temperature);
}

double FiniteElementSpecHum::est(double temperature, double salinity) const
{
    double salFactor = 1 - 5.37e-4 * salinity;
    return m_a * exp((m_b - temperature / m_d) * temperature / (temperature + m_c)) * salFactor;
}

} /* namespace Nextsim */
