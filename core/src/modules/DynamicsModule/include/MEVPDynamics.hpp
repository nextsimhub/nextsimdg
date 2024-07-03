/*!
 * @file MEVPDynamics.hpp
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

#ifndef DGCOMP
#define DGCOMP 3 // define to make red lines go away in the IDE
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

namespace Nextsim {
class MEVPDynamics : public IDynamics {
public:
    MEVPDynamics();

    std::string getName() const override { return "MEVPDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
    ModelState getStateRecursive(const OutputSpec& os) const override;
private:
    MEVPDynamicsKernel<DGCOMP> kernel;
    VPParameters params;
};
}

#endif /* DYNAMICS_HPP */
