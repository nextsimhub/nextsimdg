/*!
 * @file AdvectionDynamicsKernel.hpp
 *
 * @date Aug 1, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ADVECTIONDYNAMICSKERNEL_HPP
#define ADVECTIONDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection> class AdvectionDynamicsKernel : public CGDynamicsKernel<DGadvection> {
public:
    void update(const TimestepTime& tst) override
    {
        // Let DynamicsKernel handle the advection step
        DynamicsKernel<DGadvection, DGstressComp>::advectionAndLimits(tst);
        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressComp>::update(tst);
    }

protected:
    void updateMomentum(const TimestepTime& tst) override { }
};

} /* namespace Nextsim */

#endif /* ADVECTIONDYNAMICSKERNEL_HPP */
