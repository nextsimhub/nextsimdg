/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef LINEARFREEZING_HPP
#define LINEARFREEZING_HPP

#include "include/IFreezingPoint.hpp"
#include "include/constants.hpp"

namespace Nextsim {

//! The implementation class of the linear model of seawater freezing point.
class LinearFreezingImpl : public IFreezingPoint {
public:
    std::string getName() const override { return "LinearFreezing"; }

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
        // μ is positive, so a negative sign is needed so that the freezing point is below zero.
        return -Water::mu * sss;
    }
};

using LinearFreezing = FreezingPointImpl<LinearFreezingImpl>;
}

#endif /* LINEARFREEZING_HPP */
