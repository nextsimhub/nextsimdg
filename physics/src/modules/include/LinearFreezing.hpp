/*!
 * @file LinearFreezing.hpp
 *
 * @date Nov 10, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_LINEARFREEZING_HPP
#define SRC_INCLUDE_LINEARFREEZING_HPP

#include "IFreezingPoint.hpp"
#include "include/constants.hpp"

namespace Nextsim {

class LinearFreezing : public IFreezingPoint {
public:
    // ~LinearFreezing() = default;

    inline double operator()(double sss) const override
    {
        // μ is positive, so a negative sign is needed so that the freezing point is below zero.
        return -Water::mu * sss;
    }
};
}

#endif /* SRC_INCLUDE_LINEARFREEZING_HPP */
