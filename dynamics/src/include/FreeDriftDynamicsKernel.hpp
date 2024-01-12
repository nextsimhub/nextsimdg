/*!
 * @file FreeDriftDynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef FREEDRIFTDYNAMICSKERNEL_HPP
#define FREEDRIFTDYNAMICSKERNEL_HPP

#include "IDynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection> class FreeDriftDynamicsKernel : public IDynamicsKernel<CGdegree, DGadvection> {
using IDynamicsKernel<CGdegree, DGadvection>::NT_evp;
using IDynamicsKernel<CGdegree, DGadvection>::momentum;
using IDynamicsKernel<CGdegree, DGadvection>::hice;
using IDynamicsKernel<CGdegree, DGadvection>::cice;

protected:
    void updateMomentum(const TimestepTime& tst) override
    {
        //! Momentum
        for (size_t fdStep = 0; fdStep < NT_evp; ++fdStep) {
            momentum->freeDriftStep(FDP, NT_evp, tst.step.seconds(), hice, cice);
        }
    };

private:
    FDParameters FDP;
};
} /* namespace Nextsim */

#endif /* FREEDRIFTDYNAMICSKERNEL_HPP */
