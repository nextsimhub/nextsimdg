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

template <int DGadvection, int DGstress> class MEVPDynamicsKernel : public DynamicsKernel<DGadvection, DGstress> {
using DynamicsKernel<DGadvection, DGstress>::NT_evp;
using DynamicsKernel<DGadvection, DGstress>::momentum;
using DynamicsKernel<DGadvection, DGstress>::hice;
using DynamicsKernel<DGadvection, DGstress>::cice;

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
