/*!
 * @file VPCGDynamicsKernel.hpp
 *
 * @date 06 Dec 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef VPCGDYNAMICSKERNEL_HPP
#define VPCGDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

#include "DynamicsParameters.hpp"
#include "ParametricMap.hpp"
#include "StressUpdateStep.hpp"
#include "VPParameters.hpp"
#include "include/constants.hpp"

namespace Nextsim {

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection> class VPCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
protected:
    using DynamicsKernel<DGadvection, DGstressComp>::s11;
    using DynamicsKernel<DGadvection, DGstressComp>::s12;
    using DynamicsKernel<DGadvection, DGstressComp>::s22;
    using DynamicsKernel<DGadvection, DGstressComp>::e11;
    using DynamicsKernel<DGadvection, DGstressComp>::e12;
    using DynamicsKernel<DGadvection, DGstressComp>::e22;
    using DynamicsKernel<DGadvection, DGstressComp>::hice;
    using DynamicsKernel<DGadvection, DGstressComp>::cice;
    using DynamicsKernel<DGadvection, DGstressComp>::smesh;
    using DynamicsKernel<DGadvection, DGstressComp>::deltaT;
    using DynamicsKernel<DGadvection, DGstressComp>::stressDivergence;
    using DynamicsKernel<DGadvection, DGstressComp>::applyBoundaries;

    using CGDynamicsKernel<DGadvection>::u;
    using CGDynamicsKernel<DGadvection>::v;
    using CGDynamicsKernel<DGadvection>::xGradSeaSurfaceHeight;
    using CGDynamicsKernel<DGadvection>::yGradSeaSurfaceHeight;
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
    VPCGDynamicsKernel(StressUpdateStep<DGadvection, DGstressComp>& stressStepIn,
        const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>()
        , stressStep(stressStepIn)
        , params(reinterpret_cast<const VPParameters&>(paramsIn))
    {
    }
    virtual ~VPCGDynamicsKernel() = default;
    void update(const TimestepTime& tst) override
    {
        // Let DynamicsKernel handle the advection step
        DynamicsKernel<DGadvection, DGstressComp>::advectionAndLimits(tst);

        prepareIteration({ { hiceName, hice }, { ciceName, cice } });

        u0 = u;
        v0 = v;

        // The critical timestep for the VP solver is the advection timestep
        deltaT = tst.step.seconds();

        for (size_t mevpstep = 0; mevpstep < params.nSteps; ++mevpstep) {

            projectVelocityToStrain();

            std::array<std::reference_wrapper<DGVector<DGstressComp>>, N_TENSOR_ELEMENTS> stress
                = { s11, s12, s22 }; // Call the step function on the StressUpdateStep class
            // Call the step function on the StressUpdateStep class
            stressStep.stressUpdateHighOrder(
                params, *smesh, stress, { e11, e12, e22 }, hice, cice, deltaT);

            stressDivergence(); // Compute divergence of stress tensor

            updateMomentum(tst);

            applyBoundaries();
        }
        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressComp>::update(tst);
    }

    double getIceOceanStressElement(const std::string& name, const int i) const override
    {
        const double FOcean = params.COcean * params.rhoOcean;

        double uOcnRel = u(i) - uOcean(i);
        double vOcnRel = v(i) - vOcean(i);
        double absocn = sqrt(SQR(uOcnRel) + SQR(vOcnRel));

        if (name == uIOStressName)
            return FOcean * absocn * uOcnRel;
        else if (name == vIOStressName)
            return FOcean * absocn * vOcnRel;
        else
            return std::numeric_limits<double>::quiet_NaN();
    }

protected:
    StressUpdateStep<DGadvection, DGstressComp>& stressStep;
    const VPParameters& params;
    const double alpha = 1500.;
    const double beta = 1500.;

    // Step-initial ice velocity
    CGVector<CGdegree> u0;
    CGVector<CGdegree> v0;

    void updateMomentum(const TimestepTime& tst) override
    {

        // Update the velocity

        const double FOcean = params.COcean * params.rhoOcean;
        const double FAtm = params.CAtm * params.rhoAtm;

        //      update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
        for (int i = 0; i < u.rows(); ++i) {
            auto uOcnRel = u(i) - uOcean(i);
            auto vOcnRel = v(i) - vOcean(i);
            double absatm = sqrt(SQR(uAtmos(i)) + SQR(vAtmos(i)));
            double absocn = sqrt(
                SQR(uOcnRel) + SQR(vOcnRel)); // note that the sign of uOcnRel is irrelevant here
            auto uPrev = u(i);

            // TODO: Take the sign of lat into account for Coriolis term
            u(i) = 1.0
                / (params.rhoIce * cgH(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgA(i) * FOcean * absocn) // implicit parts
                * (params.rhoIce * cgH(i) / deltaT * (beta * u(i) + u0(i)) // pseudo-timestepping
                    + cgA(i)
                        * (FAtm * absatm * uAtmos(i) + // atm forcing
                            FOcean * absocn * uOcean(i)) // ocean forcing
                    + params.rhoIce * cgH(i) * params.fc * v(i) // Coriolis
                    - params.rhoIce * cgH(i) * PhysicalConstants::g
                        * xGradSeaSurfaceHeight(i) // sea surface
                    + dStressX(i) / pmap->lumpedcgmass(i)); // Internal stress term
            v(i) = 1.0
                / (params.rhoIce * cgH(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgA(i) * FOcean * absocn) // implicit parts
                * (params.rhoIce * cgH(i) / deltaT * (beta * v(i) + v0(i)) // pseudo-timestepping
                    + cgA(i)
                        * (FAtm * absatm * vAtmos(i) + // atm forcing
                            FOcean * absocn * vOcean(i)) // ocean forcing
                    - params.rhoIce * cgH(i) * params.fc * uPrev // Coriolis
                    - params.rhoIce * cgH(i) * PhysicalConstants::g
                        * yGradSeaSurfaceHeight(i) // sea surface
                    + dStressY(i) / pmap->lumpedcgmass(i)); // Internal stress term
        }
    }

private:
    VPCGDynamicsKernel();
};

} /* namespace Nextsim */

#endif /* VPCGDYNAMICSKERNEL_HPP */
