/*!
 * @file DummyDynamics.hpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DUMMYDYNAMICS_HPP
#define DUMMYDYNAMICS_HPP

#include "IDynamics.hpp"

#include "../../../../dynamics/src/include/DummyDynamicsKernel.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
class DummyDynamics : public IDynamics {
public:
    DummyDynamics();

    std::string getName() const override { return "DummyDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
private:
    // TODO: How to get the template parameters here?
    DummyDynamicsKernel<2, 2> kernel;
};
}

#endif /* DUMMYDYNAMICS_HPP */
