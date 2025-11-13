/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMDYNAMICS_HPP
#define BBMDYNAMICS_HPP

#include "include/BBMDynamicsKernel.hpp"
#include "include/BBMParameters.hpp"
#include "include/IDynamics.hpp"
#include "kokkos/include/KokkosBBMDynamicsKernel.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // Define to prevent errors from static analysis tools
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

namespace Nextsim {

class BBMDynamics : public IDynamics, public Configured<BBMDynamics> {
public:
    BBMDynamics();

    std::string getName() const override { return "BBMDynamics"; }
    void prepareAdvection() override;
    void update(const TimestepTime& tst) override;

    void advectField(double timestep, ModelArray& field,
        double lowerLimit = -std::numeric_limits<double>::infinity(),
        double upperLimit = std::numeric_limits<double>::infinity()) override;

#ifdef USE_KOKKOS
    void advectField(double timestep, const DeviceViewMA& field,
        double lowerLimit = -std::numeric_limits<double>::infinity(),
        double upperLimit = std::numeric_limits<double>::infinity()) override;
#endif

    void setData(const ModelState::DataMap&) override;
    void configure() override;
    ConfigMap getConfiguration() const override;

    enum {
        C_KEY,
        NU_KEY,
        YOUNG_KEY,
        P0_KEY,
        LAMBDA0_KEY,
        ALPHA_KEY,
        EXPPMAX_KEY,
        MU_KEY,
        NMAX_KEY,
        CLAB_KEY,
        NSTEPS_KEY,
        RHOI_KEY,
        RHOA_KEY,
        RHOO_KEY,
        CATM_KEY,
        COCEAN_KEY,
        FC_KEY,
        ANGLE_KEY,
    };

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap&, bool getAll);

private:
    BBMParameters params;
#ifdef USE_KOKKOS
    KokkosBBMDynamicsKernel<DGCOMP> kernel;
#else
    BBMDynamicsKernel<DGCOMP> kernel;
#endif
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICS_HPP */
