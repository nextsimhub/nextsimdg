/*!
 * @file BBMDynamicsKernel.hpp
 *
 * @date 3 Nov 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMDYNAMICSKERNEL_HPP
#define BBMDYNAMICSKERNEL_HPP

#include "IDynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection> class BBMDynamicsKernel : public IDynamicsKernel<CGdegree, DGadvection> {
using IDynamicsKernel<CGdegree, DGadvection>::NT_evp;
using IDynamicsKernel<CGdegree, DGadvection>::momentum;
using IDynamicsKernel<CGdegree, DGadvection>::hice;
using IDynamicsKernel<CGdegree, DGadvection>::cice;

protected:
    void updateMomentum(const TimestepTime& tst) override
    {
        //! Momentum
        for (size_t bbmstep = 0; bbmstep < NT_evp; ++bbmstep) {
            momentum->BBMStep(mebParams, NT_evp, tst.step.seconds(), hice, cice, damage);
        }
    };

    void setData(const std::string& name, const ModelArray& data) override
    {
        if (name == damageName) {
            DGModelArray::ma2dg(data, damage);
        } else {
            IDynamicsKernel::setData(name, data);
        }
    }
private:
    MEBParameters mebParams;
    DGVector<DGadvection> damage;


};

} /* namespace Nextsim */

#endif /* BBMDYNAMICSKERNEL_HPP */
