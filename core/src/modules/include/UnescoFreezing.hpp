/*!
 * @file UnescoFreezing.hpp
 *
 * @date Nov 10, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_UNESCOFREEZING_HPP
#define SRC_INCLUDE_UNESCOFREEZING_HPP

#include <cmath>

#include "IFreezingPoint.hpp"

namespace Nextsim {

//! The implementation class of the UNESCO model of the freezing point of
// seawater.
class UnescoFreezing : public IFreezingPoint {
    /*!
     * @brief Calculates the freezing point of seawater.
     *
     * @details Freezing point in ˚C of water with the given salinity at
     * standard pressure.
     *
     * @param sss Sea surface salinity [PSU]
     */
    inline double operator()(double sss) const override
    {
        // Fofonoff and Millard, Unesco technical papers in marine science 44, (1983)
        const double a0 = -0.0575;
        const double a1 = +1.710523e-3;
        const double a2 = -2.154996e-4;
        const double b = -7.53e-4;
        const double p0 = 0; // Zero hydrostatic pressure

        return sss * (a0 + a1 * std::sqrt(sss) + a2 * sss) + b * p0;
    }
};
}

#endif /* SRC_INCLUDE_UNESCOFREEZING_HPP */
