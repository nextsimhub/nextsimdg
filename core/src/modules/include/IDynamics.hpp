/*!
 * @file IDynamics.hpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICS_HPP
#define IDYNAMICS_HPP

#include "include/Time.hpp"

namespace Nextsim {
class IDynamics {
    IDynamics() = default;
    virtual ~IDynamics() = default;

    virtual void update(const TimestepTime& tst) = 0;
};
}

#endif /* IDYNAMICS_HPP */
