/*!
 * @file BBMDynamicsKernel.hpp
 *
 * @date 28 Apr 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef BBMDYNAMICSKERNEL_HPP
#define BBMDYNAMICSKERNEL_HPP

#include "BrittleCGDynamicsKernel.hpp"

#include "BBMStressUpdateStep.hpp"

namespace Nextsim {

template <int DGadvection> class BBMDynamicsKernel : public BrittleCGDynamicsKernel<DGadvection> {
public:
    using DynamicsKernel<DGadvection, DGstressComp>::hice;
    using DynamicsKernel<DGadvection, DGstressComp>::cice;
    using DynamicsKernel<DGadvection, DGstressComp>::smesh;
    using CGDynamicsKernel<DGadvection>::pmap;
    using BrittleCGDynamicsKernel<DGadvection>::damage;
    using BrittleCGDynamicsKernel<DGadvection>::params;
    using CGDynamicsKernel<DGadvection>::initialise;
    BBMDynamicsKernel(const BBMParameters& paramsIn)
        : BrittleCGDynamicsKernel<DGadvection>(bbmStressStep, paramsIn)
    {
    }

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
        BrittleCGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);
        bbmStressStep.setPMap(pmap.get());
        bbmStressStep.setDamage(damage);
        bbmStressStep.setCohesion(params, *smesh);
    }

private:
    // BBM stress update class
    BBMStressUpdateStep<DGadvection, DGstressComp, CGdegree> bbmStressStep;
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICSKERNEL_HPP */
