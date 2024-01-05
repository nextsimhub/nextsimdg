/*!
 * @file Dynamics.hpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP

#include "include/MEVPDynamicsKernel.hpp"
#include "include/IDynamics.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {
class Dynamics : public IDynamics {
public:
    Dynamics();

    std::string getName() const override { return "Dynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
private:
    // TODO: How to get the template parameters here?
    MEVPDynamicsKernel<2, 6> kernel;
};
}

#endif /* DYNAMICS_HPP */
