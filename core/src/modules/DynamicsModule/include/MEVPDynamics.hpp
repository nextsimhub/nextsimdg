/*!
 * @file MEVPDynamics.hpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef MEVPDYNAMICS_HPP
#define MEVPDYNAMICS_HPP

#include "include/IDamageHealing.hpp"
#include "include/IDynamics.hpp"
#include "include/MEVPDynamicsKernel.hpp"
#include "include/Module.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

extern template class Module::Module<Nextsim::IDamageHealing>;

namespace Nextsim {
class MEVPDynamics : public IDynamics, public Configured<MEVPDynamics> {
public:
    MEVPDynamics();

    std::string getName() const override { return "MEVPDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
    void configure() override;

private:
    // TODO: How to get the template parameters here?
    MEVPDynamicsKernel<6> kernel;
    VPParameters params;
};
}

#endif /* MEVPDYNAMICS_HPP */
