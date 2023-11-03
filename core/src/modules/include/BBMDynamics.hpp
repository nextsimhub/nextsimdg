/*!
 * @file BBMDynamics.hpp
 *
 * @date 3 Nov 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMDYNAMICS_HPP
#define BBMDYNAMICS_HPP

#include "IDynamics.hpp"

#include "include/BBMDynamicsKernel.hpp"

namespace Nextsim {

class BBMDynamics : public IDynamics {
public:
    BBMDynamics();

    std::string getName() const override { return "BBMDynamics"; }
    void update(const TimestepTime& tst) override;

    void setData(const ModelState::DataMap&) override;
private:
    // TODO: How to get the template parameters here?
    BBMDynamicsKernel<2, 6> kernel;
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICS_HPP */
