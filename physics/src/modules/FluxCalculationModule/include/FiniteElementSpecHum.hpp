/*
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINITEELEMENTSPECHUM_HPP
#define FINITEELEMENTSPECHUM_HPP

#include "ISpecificHumidity.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

class FiniteElementSpecHum : public ISpecificHumidity {
public:
    double operator()(double temperature, double pressure) const override;
    double operator()(double temperature, double pressure, double salinity) const override;

    std::pair<double, double> valueAndDerivative(
        double temperature, double pressure) const override;
    std::pair<double, double> valueAndDerivative(
        double temperature, double pressure, double salinity) const override;

    //! Returns a static instance already constructed to calculate specific
    //! humidity over liquid water.
    static FiniteElementSpecHum& water() { return m_water; }
    //! Returns a static instance already constructed to calculate specific
    //! humidity over ice.
    static FiniteElementSpecHum& ice() { return m_ice; }

private:
    FiniteElementSpecHum();
    // General constructor
    FiniteElementSpecHum(
        double a, double b, double c, double d, double bigA, double bigB, double bigC);

    std::pair<double, double> calculate(
        double temperature, double pressure, double salinity, bool doDeriv) const;

    double f(double temperature, double pressurePa) const;
    double est(double temperature, double salinity) const;

    const double m_a;
    const double m_b;
    const double m_c;
    const double m_d;
    const double m_bigA;
    const double m_bigB;
    const double m_bigC;
    const double m_alpha;
    const double m_beta;

    static FiniteElementSpecHum m_water;
    static FiniteElementSpecHum m_ice;

    struct Constructor {
        Constructor();
    };
    static Constructor cons;
};

// ********************************************* //
class FiniteElementSpecHum2 {
public:
    // General constructor
    constexpr FiniteElementSpecHum2(
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

    // device functions need to be defined inline
    KERNEL_IMPL_FUNCTION double operator()(double temperature, double pressure) const
    {
        return operator()(temperature, pressure, 0);
    }
    KERNEL_IMPL_FUNCTION double operator()(
        double temperature, double pressure, double salinity) const
    {
        return calculate(temperature, pressure, salinity, false).first;
    }

    KERNEL_IMPL_FUNCTION Utils::pair<double, double> valueAndDerivative(
        double temperature, double pressure) const
    {
        return valueAndDerivative(temperature, pressure, 0);
    }
    KERNEL_IMPL_FUNCTION Utils::pair<double, double> valueAndDerivative(
        double temperature, double pressure, double salinity) const
    {
        return calculate(temperature, pressure, salinity, true);
    }

private:
    KERNEL_IMPL_FUNCTION Utils::pair<double, double> calculate(
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

        return Utils::make_pair<double, double>(sphum, deriv);
    }

    // Specific humidity terms
    KERNEL_IMPL_FUNCTION double f(double temperature, double pressurePa) const
    {
        const double pressure_mb = pressurePa * 0.01;
        return 1 + m_bigA + pressure_mb * (m_bigB + m_bigC * temperature * temperature);
    }
    KERNEL_IMPL_FUNCTION double est(double temperature, double salinity) const
    {
        double salFactor = 1 - 5.37e-4 * salinity;
        return m_a * Utils::exp((m_b - temperature / m_d) * temperature / (temperature + m_c))
            * salFactor;
    }

    const double m_a;
    const double m_b;
    const double m_c;
    const double m_d;
    const double m_bigA;
    const double m_bigB;
    const double m_bigC;
    const double m_alpha;
    const double m_beta;
};

constexpr FiniteElementSpecHum2 finiteElementSpecHumWater(
    6.1121e2, 18.729, 257.87, 227.3, 7.2e-4, 3.20e-6, 5.9e-10);
constexpr FiniteElementSpecHum2 finiteElementSpecHumIce(
    6.1115e2, 23.036, 279.82, 333.7, 2.2e-4, 3.83e-6, 6.4e-10);

} /* namespace Nextsim */

#endif /* FINITEELEMENTSPECHUM_HPP */
