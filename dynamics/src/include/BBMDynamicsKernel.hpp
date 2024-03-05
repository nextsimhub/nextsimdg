/*!
 * @file BBMDynamicsKernel.hpp
 *
 * @date 3 Nov 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef BBMDYNAMICSKERNEL_HPP
#define BBMDYNAMICSKERNEL_HPP

#include "BrittleCGDynamicsKernel.hpp"

#include "BBMStressUpdateStep.hpp"

namespace Nextsim {

template <int DGadvection>
class BBMDynamicsKernel : public BrittleCGDynamicsKernel<DGadvection> {
public:
    using DynamicsKernel<CGdegree, DGadvection>::nSteps;
//using DynamicsKernel<CGdegree, DGadvection>::momentum;
    using DynamicsKernel<CGdegree, DGadvection>::hice;
    using DynamicsKernel<CGdegree, DGadvection>::cice;
    using CGDynamicsKernel<DGadvection>::pmap;
    using CGDynamicsKernel<DGadvection>::damage;
    using CGDynamicsKernel<DGadvection>::initialise;
    BBMDynamicsKernel(const DynamicsParameters& paramsIn)
        : BrittleCGDynamicsKernel<DGadvection>(bbmStressStep, paramsIn)
    {
    }

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
        CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);
        bbmStressStep.setPMap(pmap);
        bbmStressStep.setDamage(damage);
    }


private:
    //! Brittle rheology parameters
    MEBParameters mebParams;
    // BBM stress update class
    BBMStressUpdateStep<DGadvection, DGstressDegree, CGdegree> bbmStressStep;
};

} /* namespace Nextsim */

#endif /* BBMDYNAMICSKERNEL_HPP */
