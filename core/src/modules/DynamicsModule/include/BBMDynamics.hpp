/*!
 * @file BBMDynamics.hpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMDYNAMICS_HPP
#define BBMDYNAMICS_HPP

#include "include/BBMDynamicsKernel.hpp"
#include "include/IDynamics.hpp"
#include "include/MEBParameters.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // Define to prevent errors from static analysis tools
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

namespace Nextsim {

class BBMDynamics : public IDynamics {
public:
    BBMDynamics();

    std::string getName() const override { return "BBMDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

private:
    BBMDynamicsKernel<DGCOMP> kernel;
    MEBParameters params;
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICS_HPP */
