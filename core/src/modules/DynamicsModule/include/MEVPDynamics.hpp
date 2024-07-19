/*!
 * @file MEVPDynamics.hpp
 *
 * @date 18 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef MEVPDYNAMICS_HPP
#define MEVPDYNAMICS_HPP

#include "include/IDamageHealing.hpp"
#include "include/IDynamics.hpp"
#include "include/MEVPDynamicsKernel.hpp"
#include "include/Module.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // define to make red lines go away in the IDE
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

extern template class Module::Module<Nextsim::IDamageHealing>;

namespace Nextsim {
class MEVPDynamics : public IDynamics, public Configured<MEVPDynamics> {
public:
    MEVPDynamics();

    std::string getName() const override { return "MEVPDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
    ModelState getStateRecursive(const OutputSpec& os) const override;
    void configure() override;

private:
    MEVPDynamicsKernel<DGCOMP> kernel;
    VPParameters params;
};
}

#endif /* MEVPDYNAMICS_HPP */
