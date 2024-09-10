/*!
 * @file BBMDynamics.hpp
 *
 * @date Jan 5, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMDYNAMICS_HPP
#define BBMDYNAMICS_HPP

#include "include/BBMDynamicsKernel.hpp"
#include "include/IDynamics.hpp"
#include "include/MEBParameters.hpp"
#include "kokkos/include/KokkosBBMDynamicsKernel.hpp"

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
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override;

private:
#ifdef USE_KOKKOS
    KokkosBBMDynamicsKernel<DGCOMP> kernel;
#else
    BBMDynamicsKernel<DGCOMP> kernel;
#endif
    MEBParameters params;
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICS_HPP */
