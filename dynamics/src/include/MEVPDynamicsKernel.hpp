/*!
 * @file MEVPDynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef MEVPDYNAMICSKERNEL_HPP
#define MEVPDYNAMICSKERNEL_HPP

#include "VPCGDynamicsKernel.hpp"

#include "MEVPStressUpdateStep.hpp"

namespace Nextsim {

template <int DGadvection>
class MEVPDynamicsKernel : public VPCGDynamicsKernel<DGadvection> {
public:
    using CGDynamicsKernel<DGadvection>::pmap;
    using CGDynamicsKernel<DGadvection>::initialise;
    MEVPDynamicsKernel(const DynamicsParameters& paramsIn)
        : VPCGDynamicsKernel<DGadvection>(MEVPStressStep, paramsIn)
    {
    }

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
        CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);
        MEVPStressStep.setPMap(pmap);
    }

private:
    //! Rheology-Parameters
    Nextsim::VPParameters VP;
    MEVPStressUpdateStep<DGadvection, DGstressDegree, CGdegree> MEVPStressStep;
};

} /* namespace Nextsim */

#endif /* MEVPDYNAMICSKERNEL_HPP */
