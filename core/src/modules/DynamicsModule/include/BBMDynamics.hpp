/*!
 * @file BBMDynamics.hpp
 *
 * @date Jan 5, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMDYNAMICS_HPP
#define BBMDYNAMICS_HPP

#include "include/BBMDynamicsKernel.hpp"
#include "include/MEVPDynamicsKernel.hpp"
#include "include/VPParameters.hpp"
#include "include/IDynamics.hpp"

namespace Nextsim {

class BBMDynamics : public IDynamics {
public:
    BBMDynamics();

    std::string getName() const override { return "BBMDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
private:
    // TODO: How to get the template parameters here?
//    BBMDynamicsKernel<2, 6> kernel;
    // FIXME temporary use of MEVP, revert back to BBM
    MEVPDynamicsKernel<6> kernel;
    VPParameters params;

};

} /* namespace Nextsim */

#endif /* BBMDYNAMICS_HPP */
