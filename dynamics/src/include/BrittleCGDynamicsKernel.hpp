/*!
 * @file BrittleCGDynamicsKernel.hpp
 *
 * @date Feb 29, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BRITTLECGDYNAMICSKERNEL_HPP
#define BRITTLECGDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

#include "DynamicsParameters.hpp"
#include "ParametricMap.hpp"
#include "StressUpdateStep.hpp"

namespace Nextsim {

// The brittle momentum solver for CG velocity fields
template <int DGadvection>
class BrittleCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
protected:
    using DynamicsKernel<DGadvection, DGstressDegree>::nSteps;
    using DynamicsKernel<DGadvection, DGstressDegree>::s11;
    using DynamicsKernel<DGadvection, DGstressDegree>::s12;
    using DynamicsKernel<DGadvection, DGstressDegree>::s22;
    using DynamicsKernel<DGadvection, DGstressDegree>::e11;
    using DynamicsKernel<DGadvection, DGstressDegree>::e12;
    using DynamicsKernel<DGadvection, DGstressDegree>::e22;
    using DynamicsKernel<DGadvection, DGstressDegree>::hice;
    using DynamicsKernel<DGadvection, DGstressDegree>::cice;
    using DynamicsKernel<DGadvection, DGstressDegree>::smesh;
    using DynamicsKernel<DGadvection, DGstressDegree>::deltaT;
    using DynamicsKernel<DGadvection, DGstressDegree>::stressDivergence;
    using DynamicsKernel<DGadvection, DGstressDegree>::applyBoundaries;

    using CGDynamicsKernel<DGadvection>::u;
    using CGDynamicsKernel<DGadvection>::v;
    using CGDynamicsKernel<DGadvection>::uAtmos;
    using CGDynamicsKernel<DGadvection>::vAtmos;
    using CGDynamicsKernel<DGadvection>::uOcean;
    using CGDynamicsKernel<DGadvection>::vOcean;
    using CGDynamicsKernel<DGadvection>::prepareIteration;
    using CGDynamicsKernel<DGadvection>::projectVelocityToStrain;
    using CGDynamicsKernel<DGadvection>::cgH;
    using CGDynamicsKernel<DGadvection>::cgA;
    using CGDynamicsKernel<DGadvection>::dStressX;
    using CGDynamicsKernel<DGadvection>::dStressY;
    using CGDynamicsKernel<DGadvection>::pmap;
public:
    BrittleCGDynamicsKernel(StressUpdateStep<DGadvection, DGstressDegree>& stressStepIn, const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>(),
          stressStep(stressStepIn),
          params(reinterpret_cast<const VPParameters&>(paramsIn))
    {
    }
    virtual ~BrittleCGDynamicsKernel() = default;
    void update(const TimestepTime& tst) override {
        // The timestep for the brittle solver is the solver subtimestep
        deltaT = tst.step.seconds() / nSteps;

        for (size_t mevpstep = 0; mevpstep < nSteps; ++mevpstep) {

            projectVelocityToStrain();

            std::array<DGVector<DGstressDegree>, 3> stress = { s11, s12, s22 };
            // Call the step function on the StressUpdateStep class
            stressStep.stressUpdateHighOrder(params,
                    *smesh,
                    stress, { e11, e12, e22 },
                    hice, cice,
                    deltaT);

            s11 = stress[I11];
            s12 = stress[I12];
            s22 = stress[I22];

            dStressX = 0;
            dStressY = 0;

            stressDivergence(1.0);

            // Common brittle parts of the momentum solver.
        }
    }

};
}

#endif /* BRITTLECGDYNAMICSKERNEL_HPP */
