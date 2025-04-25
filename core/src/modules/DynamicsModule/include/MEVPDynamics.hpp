/*!
 * @file MEVPDynamics.hpp
 *
 * @date 27 Mar 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef MEVPDYNAMICS_HPP
#define MEVPDYNAMICS_HPP

#include "include/IDamageHealing.hpp"
#include "include/IDynamics.hpp"
#include "include/MEVPDynamicsKernel.hpp"
#include "include/NextsimModule.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // Define to prevent errors from static analysis tools
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
    };

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap&, bool getAll);

private:
    VPParameters params;
    MEVPDynamicsKernel<DGCOMP> kernel;
};
}

#endif /* MEVPDYNAMICS_HPP */
