/*!
 * @file FreeDriftDynamics.hpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FREEDRIFTDYNAMICS_HPP
#define FREEDRIFTDYNAMICS_HPP

#include "include/FreeDriftDynamicsKernel.hpp"
#include "include/IDynamics.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
class FreeDriftDynamics : public IDynamics {
public:
    FreeDriftDynamics()
        : IDynamics()
    {
    }

    std::string getName() const override { return "MEVPDynamics"; }
    void update(const TimestepTime& tst) override { kernel.update(tst); }

private:
    // TODO: How to get the template parameters here?
    FreeDriftDynamicsKernel<6> kernel;
};
}

#endif /* FREEDRIFTDYNAMICS_HPP */
