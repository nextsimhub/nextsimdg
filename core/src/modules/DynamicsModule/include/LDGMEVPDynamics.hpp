/*!
 * @file LDGMEVPDynamics.hpp
 *
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef LDGMEVPDYNAMICS_HPP
#define LDGMEVPDYNAMICS_HPP

#include "include/IDamageHealing.hpp"
#include "include/IDynamics.hpp"
#include "include/LDGMEVPDynamicsKernel.hpp"
#include "include/NextsimModule.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // Define to prevent errors from static analysis tools
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

#ifndef LDGCOMP
#define LDGCOMP 4 // Define to prevent errors from static analysis tools
#error "Number of DG components of the LDG velocity (LDGCOMP) not defined"
#endif

extern template class Module::Module<Nextsim::IDamageHealing>;

namespace Nextsim {
/*!
 * mEVP dynamics with the LDG discretization of the momentum equation. Cartesian grids only,
 * no periodic boundaries, no MPI and no Kokkos.
 */
class LDGMEVPDynamics : public IDynamics, public Configured<LDGMEVPDynamics> {
public:
    LDGMEVPDynamics();

    std::string getName() const override { return "LDGMEVPDynamics"; }
    void prepareAdvection() override;
    void update(const TimestepTime& tst) override;

    void advectField(FloatType timestep, ModelArray& field,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity()) override;

#ifdef USE_KOKKOS
    void advectField(FloatType timestep, const DeviceViewMA& field,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity()) override;
#endif

    void setData(const ModelState::DataMap&) override;
    void configure() override;
    ConfigMap getConfiguration() const override;

    enum {
        PSTAR_KEY,
        DELTA_KEY,
        C_KEY,
        NSTEPS_KEY,
        RHOI_KEY,
        RHOA_KEY,
        RHOO_KEY,
        CATM_KEY,
        COCEAN_KEY,
        FC_KEY,
        ANGLE_KEY,
        REPLACEMENT_KEY,
        ALPHA_KEY,
        BETA_KEY,
        LDGA_KEY,
        LDGB_KEY,
        PCGTOL_KEY,
        PCGMAXIT_KEY,
    };

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap&, bool getAll);

private:
    VPParameters params;
    LDGParameters ldgParams;
    LDGMEVPDynamicsKernel<DGCOMP, LDGCOMP> kernel;
};
}

#endif /* LDGMEVPDYNAMICS_HPP */
