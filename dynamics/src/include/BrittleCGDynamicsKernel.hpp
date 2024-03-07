/*!
 * @file BrittleCGDynamicsKernel.hpp
 *
 * @date Feb 29, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BRITTLECGDYNAMICSKERNEL_HPP
#define BRITTLECGDYNAMICSKERNEL_HPP

#include "CGDynamicsKernel.hpp"

#include "MEBParameters.hpp"
#include "ParametricMap.hpp"
#include "StressUpdateStep.hpp"

namespace Nextsim {

// Degrees to radians as a hex float
static const double radians = 0x1.1df46a2529d39p-6;

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
    using CGDynamicsKernel<DGadvection>::dirichletZero;
    using CGDynamicsKernel<DGadvection>::cgH;
    using CGDynamicsKernel<DGadvection>::cgA;
    using CGDynamicsKernel<DGadvection>::dStressX;
    using CGDynamicsKernel<DGadvection>::dStressY;
    using CGDynamicsKernel<DGadvection>::pmap;
public:
    BrittleCGDynamicsKernel(StressUpdateStep<DGadvection, DGstressDegree>& stressStepIn, const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>(),
          stressStep(stressStepIn),
          params(reinterpret_cast<const MEBParameters&>(paramsIn))
    {
    }
    virtual ~BrittleCGDynamicsKernel() = default;
    void update(const TimestepTime& tst) override {
        // The timestep for the brittle solver is the solver subtimestep
        deltaT = tst.step.seconds() / nSteps;

        for (size_t subStep = 0; subStep < nSteps; ++subStep) {

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

            dStressX.zero();
            dStressY.zero();

            stressDivergence(1.0);

            updateMomentum(tst);

            applyBoundaries();

            // Land mask
        }
        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressDegree>::update(tst);

    }

    void setData(const std::string& name, const ModelArray& data) override
    {
        if (name == damageName) {
            DGModelArray::ma2dg(data, damage);
        } else {
            CGDynamicsKernel<DGadvection>::setData(name, data);
        }
    }

protected:
    CGVector<CGdegree> avgX;
    CGVector<CGdegree> avgY;

    StressUpdateStep<DGadvection, DGstressDegree>& stressStep;
    const MEBParameters& params;

    DGVector<DGadvection> damage;

    // Common brittle parts of the momentum solver.
    void updateMomentum(const TimestepTime& tst) override
    {
        static const double cosOceanAngle = cos(radians * params.ocean_turning_angle);
        static const double sinOceanAngle = sin(radians * params.ocean_turning_angle);

#pragma omp parallel for
        for (size_t i = 0; i < u.rows(); ++i) {
            // FIXME dte_over_mass should include snow in the total mass
            const double dteOverMass = deltaT / (params.rho_ice * cgH(i));
            // Memoized initial velocity values
            const double uIce = u(i);
            const double vIce = v(i);

            const double cPrime = cgA(i) * params.F_ocean * std::hypot(uOcean(i) - uIce, vOcean(i) - vIce);

            // FIXME grounding term tauB = cBu[i] / std::hypot(uIce, vIce) + u0
            const double tauB = 0.;
            const double alpha = 1 + dteOverMass * (cPrime * cosOceanAngle + tauB);
            /* FIXME latitude needed for spherical cases
             * const double beta = deltaT * params.fc +
             * dteOverMass * cPrime * std::copysign(sinOceanAngle, lat[i]);
             */
            const double beta = deltaT * params.fc + dteOverMass * cPrime * sinOceanAngle;
            const double rDenom = 1 / (SQR(alpha) + SQR(beta));

            // Atmospheric drag
            const double dragAtm = cgA(i) * params.F_atm * std::hypot(uAtmos(i), vAtmos(i));
            const double tauX = dragAtm * uAtmos(i) +
                    cPrime * (uOcean(i) * cosOceanAngle - vOcean(i) * sinOceanAngle);
            const double tauY = dragAtm * vAtmos(i) +
                    cPrime * (vOcean(i) * cosOceanAngle + uOcean(i) * sinOceanAngle);

            // Stress gradient
            const double gradX = dStressX(i) / pmap->lumpedcgmass(i);
            const double gradY = dStressY(i) / pmap->lumpedcgmass(i);

            u(i) = alpha * uIce + beta * vIce
                    + dteOverMass * (alpha * (gradX + tauX) + beta * (gradY + tauY));
            u(i) *= rDenom;

            v(i) = alpha * vIce - beta *uIce
                    + dteOverMass * (alpha * (gradY + tauY) + beta * (gradX + tauX));
            v(i) *= rDenom;
        }
        dirichletZero(u);
        dirichletZero(v);

        // Mask the land on the CG grid, using the finite volume landmask
        static const size_t cgRowLength = CGdegree * smesh->nx + 1;
#pragma omp parallel for
        for (size_t eid = 0; eid < smesh->nelements; ++eid) {
            if (smesh->landmask[eid] == 0) {
                const size_t ex = eid % smesh->nx;
                const size_t ey = eid % smesh->ny;
                // Loop over CG elements for this finite volume grid cell
                for (size_t jy = 0; jy <= CGdegree; ++jy) {
                    for (size_t jx = 0; jx <= CGdegree; ++ jx) {
                        const size_t cgi = cgRowLength * (CGdegree * ey + jy) + CGdegree * ex + jx;
                        u(cgi) = 0.;
                        v(cgi) = 0.;
                    }
                }
            }
        }

        std::cout << __FILE__ << " done" << std::endl;

        // Calculate the contribution to the average velocity
        avgX += u / nSteps;
        avgY += v / nSteps;
    }

};
}

#endif /* BRITTLECGDYNAMICSKERNEL_HPP */
