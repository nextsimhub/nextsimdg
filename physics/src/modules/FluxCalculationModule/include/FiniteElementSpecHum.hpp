/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie<robert.jendersie@ovgu.de>
 */

#ifndef FINITEELEMENTSPECHUM_HPP
#define FINITEELEMENTSPECHUM_HPP

#include "ISpecificHumidity.hpp"
#include "include/KernelAlternatives.hpp"

namespace Nextsim {

class FiniteElementSpecHum {
public:
    // device functions need to be defined inline
    KERNEL_IMPL_FUNCTION FloatType operator()(FloatType temperature, FloatType pressure) const
    {
        return operator()(temperature, pressure, 0);
    }
    KERNEL_IMPL_FUNCTION FloatType operator()(
        FloatType temperature, FloatType pressure, FloatType salinity) const
    {
        return calculate<false>(temperature, pressure, salinity);
    }

    KERNEL_IMPL_FUNCTION Utils::pair<FloatType, FloatType> valueAndDerivative(
        FloatType temperature, FloatType pressure) const
    {
        return valueAndDerivative(temperature, pressure, 0);
    }
    KERNEL_IMPL_FUNCTION Utils::pair<FloatType, FloatType> valueAndDerivative(
        FloatType temperature, FloatType pressure, FloatType salinity) const
    {
        return calculate<true>(temperature, pressure, salinity);
    }

    //! Returns a static instance already constructed to calculate specific
    //! humidity over liquid water.
    static FiniteElementSpecHum& water() { return m_water; }
    //! Returns a static instance already constructed to calculate specific
    //! humidity over ice.
    static FiniteElementSpecHum& ice() { return m_ice; }

private:
    FiniteElementSpecHum(FloatType a, FloatType b, FloatType c, FloatType d, FloatType bigA,
        FloatType bigB, FloatType bigC)
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

    template <bool doDeriv>
    KERNEL_IMPL_FUNCTION auto calculate(
        FloatType temperature, FloatType pressure, FloatType salinity) const
        -> std::conditional_t<doDeriv, Utils::pair<FloatType, FloatType>, FloatType>
    {
        const FloatType estCalc = est(temperature, salinity);
        const FloatType fCalc = f(temperature, pressure);
        const FloatType sphum = m_alpha * fCalc * estCalc / (pressure - m_beta * fCalc * estCalc);

        if constexpr (doDeriv) {
            const FloatType df_dT = 2 * m_bigC * m_bigB * temperature;
            FloatType numerator = m_b * m_c * m_d - temperature * (2 * m_c + temperature);
            FloatType sqrtDenom = m_c + temperature;
            FloatType denominator = m_d * sqrtDenom * sqrtDenom;
            const FloatType dest_dT = numerator / denominator * estCalc;
            numerator = m_alpha * pressure * (fCalc * dest_dT + estCalc * df_dT);
            sqrtDenom = pressure - m_beta * estCalc * fCalc;
            denominator = sqrtDenom * sqrtDenom;

            const FloatType deriv = numerator / denominator;

            return Utils::pair<FloatType, FloatType>(sphum, deriv);
        } else {
            return sphum;
        }
    }

    // Specific humidity terms
    KERNEL_IMPL_FUNCTION FloatType f(FloatType temperature, FloatType pressurePa) const
    {
        const FloatType pressure_mb = pressurePa * 0.01;
        return 1 + m_bigA + pressure_mb * (m_bigB + m_bigC * temperature * temperature);
    }
    KERNEL_IMPL_FUNCTION FloatType est(FloatType temperature, FloatType salinity) const
    {
        FloatType salFactor = 1 - 5.37e-4 * salinity;
        return m_a * Utils::exp((m_b - temperature / m_d) * temperature / (temperature + m_c))
            * salFactor;
    }

    const FloatType m_a;
    const FloatType m_b;
    const FloatType m_c;
    const FloatType m_d;
    const FloatType m_bigA;
    const FloatType m_bigB;
    const FloatType m_bigC;
    const FloatType m_alpha;
    const FloatType m_beta;

    static FiniteElementSpecHum m_water;
    static FiniteElementSpecHum m_ice;
};

} /* namespace Nextsim */

#endif /* FINITEELEMENTSPECHUM_HPP */
