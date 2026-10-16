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
class UnescoFreezingImpl {
public:
    std::string getName() const { return "UnescoFreezing"; }

    /*!
     * @brief Calculates the freezing point of seawater.
     *
     * @details Freezing point in ˚C of water with the given salinity at
     * standard pressure.
     *
     * @param sss Sea surface salinity [PSU]
     */
    KERNEL_IMPL_FUNCTION FloatType calculate(FloatType sss) const
    {
        // Fofonoff and Millard, Unesco technical papers in marine science 44, (1983)
        constexpr FloatType a0 = -0.0575;
        constexpr FloatType a1 = +1.710523e-3;
        constexpr FloatType a2 = -2.154996e-4;
        constexpr FloatType b = -7.53e-4;
        constexpr FloatType p0 = 0; // Zero hydrostatic pressure

        return sss * (a0 + a1 * Utils::sqrt(sss) + a2 * sss) + b * p0;
    }
};

using UnescoFreezing = FreezingPointImpl<UnescoFreezingImpl>;
}

#endif /* UNESCOFREEZING_HPP */
