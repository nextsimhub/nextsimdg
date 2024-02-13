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
    using DynamicsKernel<DGadvection, DGstressDegree>::NT_evp;
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
    using CGDynamicsKernel<DGadvection>::projectVelocityToStrain;
    using CGDynamicsKernel<DGadvection>::cgH;
    using CGDynamicsKernel<DGadvection>::cgA;
    using CGDynamicsKernel<DGadvection>::dStressX;
    using CGDynamicsKernel<DGadvection>::dStressY;
    using CGDynamicsKernel<DGadvection>::pmap;
    VPCGDynamicsKernel(StressUpdateStep<DGadvection, DGstressDegree>& stressStepIn, const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>(),
          stressStep(stressStepIn),
          params(reinterpret_cast<const VPParameters&>(paramsIn))
    {
    }
    virtual ~VPCGDynamicsKernel() = default;
    void update(const TimestepTime& tst) override
    {
        // Let DynamicsKernel handle the advection step
        DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits(tst);

        u0 = u;
        v0 = v;

        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {

            projectVelocityToStrain();

            // FIXME Why can I not use the typedef from StressUpdateStep?
            /*StressUpdateStep<DGadvection, DGstressDegree>::SymmetricTensorVector*/
            std::array<DGVector<DGstressDegree>, 3> stress = { s11, s12, s22 };
            // Call the step function on the StressUpdateStep class
            stressStep.stressUpdateHighOrder(params,
                    *smesh,
                    stress, { e11, e12, e22 },
                    hice, cice,
                    deltaT);

            double stressScale = 1.0; // 2nd-order Stress term has different scaling with the EarthRadius
            if (smesh->CoordinateSystem == Nextsim::SPHERICAL)
                stressScale = 1.0 / Nextsim::EarthRadius / Nextsim::EarthRadius;

            stressDivergence(stressScale); // Compute divergence of stress tensor

            updateMomentum(tst);

            applyBoundaries();
        }
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
        double SC = 1.0;///(1.0-pow(1.0+1.0/beta,-1.0*NT_evp));

        //      update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
        for (int i = 0; i < u.rows(); ++i) {
            auto uOcnRel = uOcean(i) - u(i); // note the reversed sign compared to the v component
            auto vOcnRel = v(i) - vOcean(i);
            double absatm = sqrt(SQR(uAtmos(i)) + SQR(vAtmos(i)));
            double absocn = sqrt(SQR(uOcnRel) + SQR(vOcnRel)); // note that the sign of uOcnRel is irrelevant here

            u(i) = (1.0
                    / (params.rho_ice * cgH(i) / deltaT * (1.0 + beta) // implicit parts
                            + cgA(i) * params.F_ocean
                            * absocn ) // implicit parts
                            * (params.rho_ice * cgH(i) / deltaT * (beta * u(i) + u0(i)) // pseudo-timestepping
                                    + cgA(i)
                                    * (params.F_atm * absatm * uAtmos(i) + // atm forcing
                                            params.F_ocean * absocn * SC
                                            * uOcean(i)) // ocean forcing
                                            + params.rho_ice * cgH(i) * params.fc
                                            * vOcnRel // cor + surface
                                            + dStressX(i)/pmap->lumpedcgmass(i)
                            ));
            v(i) = (1.0
                    / (params.rho_ice * cgH(i) / deltaT * (1.0 + beta) // implicit parts
                            + cgA(i) * params.F_ocean
                            * absocn ) // implicit parts
                            * (params.rho_ice * cgH(i) / deltaT * (beta * v(i) + v0(i)) // pseudo-timestepping
                                    + cgA(i)
                                    * (params.F_atm * absatm * vAtmos(i) + // atm forcing
                                            params.F_ocean * absocn * SC
                                            * vOcean(i)) // ocean forcing
                                            + params.rho_ice * cgH(i) * params.fc
                                            * uOcnRel // here the reversed sign of uOcnRel is used
                                            + dStressY(i)/pmap->lumpedcgmass(i)
                            ));
        }

    }
private:
    VPCGDynamicsKernel();
};

} /* namespace Nextsim */

#endif /* VPCGDYNAMICSKERNEL_HPP */
