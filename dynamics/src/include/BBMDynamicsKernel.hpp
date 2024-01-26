/*!
 * @file BBMDynamicsKernel.hpp
 *
 * @date 3 Nov 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef BBMDYNAMICSKERNEL_HPP
#define BBMDYNAMICSKERNEL_HPP

#include "DynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection> class BBMDynamicsKernel : public DynamicsKernel<CGdegree, DGadvection> {
using DynamicsKernel<CGdegree, DGadvection>::NT_evp;
using DynamicsKernel<CGdegree, DGadvection>::momentum;
using DynamicsKernel<CGdegree, DGadvection>::hice;
using DynamicsKernel<CGdegree, DGadvection>::cice;

public:
    void setData(const std::string& name, const ModelArray& data) override
    {
        if (name == damageName) {
            DGModelArray::ma2dg(data, damage);
        } else {
            DynamicsKernel<CGdegree, DGadvection>::setData(name, data);
        }
    }

protected:
    void updateMomentum(const TimestepTime& tst) override
    {
        //! Momentum
        for (size_t bbmstep = 0; bbmstep < NT_evp; ++bbmstep) {
            momentum->BBMStep(mebParams, NT_evp, tst.step.seconds(), hice, cice, damage);
        }
    };

private:
    //! Brittle rheology parameters
    MEBParameters mebParams;
    DGVector<DGadvection> damage;
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICSKERNEL_HPP */
