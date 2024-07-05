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
#include "VPParameters.hpp"
#include "stopwatch.hpp"

namespace Nextsim {

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection> class VPCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
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
    VPCGDynamicsKernel(StressUpdateStep<DGadvection, DGstressDegree>& stressStepIn,
        const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>()
        , stressStep(stressStepIn)
        , params(static_cast<const VPParameters&>(paramsIn))
    {
    }
    virtual ~VPCGDynamicsKernel() = default;
    void update(const TimestepTime& tst) override
    {
        static PerfTimerAlt<true> timerMevp("mevp");
        static PerfTimerAlt<MEASURE_DETAILED> timerProj("proj");
        static PerfTimerAlt<MEASURE_DETAILED> timerStress("stress");
        static PerfTimerAlt<MEASURE_DETAILED> timerDivergence("div");
        static PerfTimerAlt<MEASURE_DETAILED> timerMomentum("momentum");
        static PerfTimerAlt<MEASURE_DETAILED> timerBoundary("bc");

        // Let DynamicsKernel handle the advection step
        DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits(tst);

        prepareIteration({ { hiceName, hice }, { ciceName, cice } });

        timerMevp.start();
        u0 = u;
        v0 = v;

        // The critical timestep for the VP solver is the advection timestep
        deltaT = tst.step.seconds();
        
        for (size_t mevpstep = 0; mevpstep < nSteps; ++mevpstep) {
            timerProj.start();
            projectVelocityToStrain();
            timerProj.stop();

            timerStress.start();
            std::array<std::reference_wrapper<DGVector<DGstressDegree>>, N_TENSOR_ELEMENTS> stress
                = { s11, s12, s22 }; // Call the step function on the StressUpdateStep class
            // Call the step function on the StressUpdateStep class
            stressStep.stressUpdateHighOrder(
                params, *smesh, stress, { e11, e12, e22 }, hice, cice, deltaT);
            timerStress.stop();

            timerDivergence.start();
            stressDivergence(); // Compute divergence of stress tensor
            timerDivergence.stop();

            timerMomentum.start();
            updateMomentum(tst);
            timerMomentum.stop();

            timerBoundary.start();
            applyBoundaries();
            timerBoundary.stop();
        }
        timerMevp.stop();

        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressDegree>::update(tst);
    }

protected:
    StressUpdateStep<DGadvection, DGstressDegree>& stressStep;
    const VPParameters& params;
    const double alpha = 1500.;
    const double beta = 1500.;

    // Step-initial ice velocity
    CGVector<CGdegree> u0;
    CGVector<CGdegree> v0;

    void updateMomentum(const TimestepTime& tst) override
    {

        // Update the velocity
        double SC = 1.0; ///(1.0-pow(1.0+1.0/beta,-1.0*nSteps));

        //      update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
        for (int i = 0; i < u.rows(); ++i) {
            auto uOcnRel = uOcean(i) - u(i); // note the reversed sign compared to the v component
            auto vOcnRel = v(i) - vOcean(i);
            double absatm = sqrt(SQR(uAtmos(i)) + SQR(vAtmos(i)));
            double absocn = sqrt(
                SQR(uOcnRel) + SQR(vOcnRel)); // note that the sign of uOcnRel is irrelevant here

            u(i) = (1.0
                / (params.rho_ice * cgH(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgA(i) * params.F_ocean * absocn) // implicit parts
                * (params.rho_ice * cgH(i) / deltaT * (beta * u(i) + u0(i)) // pseudo-timestepping
                    + cgA(i)
                        * (params.F_atm * absatm * uAtmos(i) + // atm forcing
                            params.F_ocean * absocn * SC * uOcean(i)) // ocean forcing
                    + params.rho_ice * cgH(i) * params.fc * vOcnRel // cor + surface
                    + dStressX(i) / pmap->lumpedcgmass(i)));
            v(i) = (1.0
                / (params.rho_ice * cgH(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgA(i) * params.F_ocean * absocn) // implicit parts
                * (params.rho_ice * cgH(i) / deltaT * (beta * v(i) + v0(i)) // pseudo-timestepping
                    + cgA(i)
                        * (params.F_atm * absatm * vAtmos(i) + // atm forcing
                            params.F_ocean * absocn * SC * vOcean(i)) // ocean forcing
                    + params.rho_ice * cgH(i) * params.fc
                        * uOcnRel // here the reversed sign of uOcnRel is used
                    + dStressY(i) / pmap->lumpedcgmass(i)));
        }
    }

private:
    VPCGDynamicsKernel();
};

} /* namespace Nextsim */

#endif /* VPCGDYNAMICSKERNEL_HPP */
