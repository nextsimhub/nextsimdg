/*!
 * @file DummyDynamics.hpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYDYNAMICS_HPP
#define DUMMYDYNAMICS_HPP

#include "include/IDynamics.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
class DummyDynamics : public IDynamics {
public:
    DummyDynamics() = default;

    std::string getName() const override { return "DummyDynamics"; }
    void update(const TimestepTime& tst) override { };

};
}

#endif /* DUMMYDYNAMICS_HPP */
