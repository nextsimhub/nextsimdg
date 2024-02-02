/*!
 * @file VPCGDynamicsKernel.hpp
 *
 * @date Feb 2, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef VPCGDYNAMICSKERNEL_HPP
#define VPCGDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

#include "DynamicsParameters.hpp"
#include "ParametricMap.hpp"
#include "StressUpdateStep.hpp"

namespace Nextsim {

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection>
class VPCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
public:
    VPCGDynamicsKernel(ParametricMap<CGdegree>* pmapIn, StressUpdateStep<DGadvection, DGstressDegree>& stressStepIn, const DynamicsParameters& paramsIn)
        : pmap(pmapIn),
          stressStep(stressStepIn),
          params(paramsIn)
    {
    }
    virtual ~VPCGDynamicsKernel() = default;
    void update(const TimestepTime& tst) override
    {
        // Let DynamicsKernel handle the advection step
        DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits(tst);

        // Call the step function on the StressUpdateStep class
        stressStep.stressUpdateHighOrder(params,
                *smesh,
                { s11, s12, s22 }, { e11, e12, e22 },
                hice, cice,
                deltaT);

        updateMomentum(tst);

        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressDegree>::update(tst);

    }
protected:
    ParametricMap<CGdegree>* pmap;
    StressUpdateStep<DGadvection, DGstressDegree>& stressStep;
    const DynamicsParameters& params;
private:
    VPCGDynamicsKernel();
};

} /* namespace Nextsim */

#endif /* VPCGDYNAMICSKERNEL_HPP */
