/*!
 * @file DummyDynamics.hpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYDYNAMICS_HPP
#define DUMMYDYNAMICS_HPP

#include "IDynamics.hpp"

namespace Nextsim {
class DummyDynamics : public IDynamics{
    void update(const TimestepTime& tst) override { }
};
}

#endif /* DUMMYDYNAMICS_HPP */
