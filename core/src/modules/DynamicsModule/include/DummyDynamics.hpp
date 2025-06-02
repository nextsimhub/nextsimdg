/*!
 * @file DummyDynamics.hpp
 *
 * @date 02 Jun 2025
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
    void update(const TimestepTime& tst) override {};
    void setData(const ModelState::DataMap& state) override
    {
        IDynamics::setData(state);

        // Set the ice velocity to zero so that we don't trip the PrognosticData field checks.
        uice = 0.;
        vice = 0.;
    }
};
}

#endif /* DUMMYDYNAMICS_HPP */
