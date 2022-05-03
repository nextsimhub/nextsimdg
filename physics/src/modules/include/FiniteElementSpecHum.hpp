/*
 * @file FiniteElementSpecHum.hpp
 *
 * @date May 3, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINITEELEMENTSPECHUM_HPP
#define FINITEELEMENTSPECHUM_HPP

#include "ISpecificHumidity.hpp"

#include "include/Configured.hpp"

namespace Nextsim {

class FiniteElementSpecHum: public ISpecificHumidity, public Configured<FiniteElementSpecHum> {
public:

    void configure() override;

    double operator()(double temperature, double pressure) const override;
    double operator()(double temperature, double pressure, double salinity) const override;

    std::pair<double, double> valueAndDerivative(double temperature, double pressure) const override;
    std::pair<double, double> valueAndDerivative(double temperature,
            double pressure, double salinity) const override;

private:
    FiniteElementSpecHum();
    // General constructor
    FiniteElementSpecHum(double a, double b, double c, double d, double bigA, double bigB, double bigC);

    std::pair<double, double> calculate(double temperature, double pressure, double salinity, bool doDeriv) const;

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

};

} /* namespace Nextsim */

#endif /* FINITEELEMENTSPECHUM_HPP */
