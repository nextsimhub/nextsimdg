/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef UNESCOFREEZING_HPP
#define UNESCOFREEZING_HPP

#include <cmath>

#include "include/IFreezingPoint.hpp"

namespace Nextsim {

//! The implementation class of the UNESCO model of the freezing point of
// seawater.
class UnescoFreezingImpl : public IFreezingPoint {
public:
    std::string getName() const override { return "UnescoFreezing"; }

    /*!
     * @brief Calculates the freezing point of seawater.
     *
     * @details Freezing point in ˚C of water with the given salinity at
     * standard pressure.
     *
     * @param sss Sea surface salinity [PSU]
     */
    KERNEL_IMPL_FUNCTION double calculate(double sss) const
    {
        // Fofonoff and Millard, Unesco technical papers in marine science 44, (1983)
        constexpr double a0 = -0.0575;
        constexpr double a1 = +1.710523e-3;
        constexpr double a2 = -2.154996e-4;
        constexpr double b = -7.53e-4;
        constexpr double p0 = 0; // Zero hydrostatic pressure

        return sss * (a0 + a1 * Utils::sqrt(sss) + a2 * sss) + b * p0;
    }
};

using UnescoFreezing = FreezingPointImpl<UnescoFreezingImpl>;
}

#endif /* UNESCOFREEZING_HPP */
