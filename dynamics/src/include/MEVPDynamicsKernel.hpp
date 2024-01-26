/*!
 * @file MEVPDynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef MEVPDYNAMICSKERNEL_HPP
#define MEVPDYNAMICSKERNEL_HPP

#include "DynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection> class MEVPDynamicsKernel : public DynamicsKernel<CGdegree, DGadvection> {
using DynamicsKernel<CGdegree, DGadvection>::NT_evp;
using DynamicsKernel<CGdegree, DGadvection>::momentum;
using DynamicsKernel<CGdegree, DGadvection>::hice;
using DynamicsKernel<CGdegree, DGadvection>::cice;

protected:
    void updateMomentum(const TimestepTime& tst) override
    {
        //! Momentum
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {
            momentum->mEVPStep(VP, NT_evp, alpha, beta, tst.step.seconds(), hice, cice);
        }
    };

private:
    //! Rheology-Parameters
    Nextsim::VPParameters VP;
    //! MEVP parameters
    double alpha = 1500.0;
    double beta = 1500.0;
};

} /* namespace Nextsim */

#endif /* MEVPDYNAMICSKERNEL_HPP */
