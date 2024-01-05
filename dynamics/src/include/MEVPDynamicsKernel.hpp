/*!
 * @file DynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef MEVPDYNAMICSKERNEL_HPP
#define MEVPDYNAMICSKERNEL_HPP

#include "IDynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection> class MEVPDynamicsKernel : public IDynamicsKernel<CGdegree, DGadvection> {
using IDynamicsKernel<CGdegree, DGadvection>::NT_evp;
using IDynamicsKernel<CGdegree, DGadvection>::momentum;
using IDynamicsKernel<CGdegree, DGadvection>::hice;
using IDynamicsKernel<CGdegree, DGadvection>::cice;

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
