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
template <int DGadvection> class BrittleCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
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
    using DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits;
    using DynamicsKernel<DGadvection, DGstressDegree>::dgtransport;

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

    double cosOceanAngle, sinOceanAngle;

public:
    BrittleCGDynamicsKernel(StressUpdateStep<DGadvection, DGstressDegree>& stressStepIn,
        const DynamicsParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>()
        , stressStep(stressStepIn)
        , params(reinterpret_cast<const MEBParameters&>(paramsIn))
        , stresstransport(nullptr)
    {
    }
    virtual ~BrittleCGDynamicsKernel() = default;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
        CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

        //! Initialize stress transport
        stresstransport = new Nextsim::DGTransport<DGstressDegree>(*smesh);
        stresstransport->settimesteppingscheme("rk2");

        damage.resize_by_mesh(*smesh);
        avgU.resize_by_mesh(*smesh);
        avgV.resize_by_mesh(*smesh);

        cosOceanAngle = cos(radians * params.ocean_turning_angle);
        sinOceanAngle = sin(radians * params.ocean_turning_angle);
    }

    // The brittle rheologies use avgU and avgV to do the advection, not u and v, like mEVP
    void prepareAdvection() override
    {
        dgtransport->prepareAdvection(avgU, avgV);
    }

    void update(const TimestepTime& tst) override
    {

        // Let DynamicsKernel handle the advection step
        advectionAndLimits(tst);

        //! Perform transport step for stress
        stresstransport->prepareAdvection(avgU, avgV);
        stresstransport->step(tst.step.seconds(), s11);
        stresstransport->step(tst.step.seconds(), s12);
        stresstransport->step(tst.step.seconds(), s22);

        // Transport and limits for damage
        dgtransport->step(tst.step.seconds(), damage);
        Nextsim::LimitMax(damage, 1.0);
        Nextsim::LimitMin(damage, 1e-12);

        prepareIteration({ { hiceName, hice }, { ciceName, cice } });

        // The timestep for the brittle solver is the solver subtimestep
        deltaT = tst.step.seconds() / nSteps;

        avgU.zero();
        avgV.zero();

        for (size_t subStep = 0; subStep < nSteps; ++subStep) {

            projectVelocityToStrain();

            std::array<std::reference_wrapper<DGVector<DGstressDegree>>, N_TENSOR_ELEMENTS> stress
                = { s11, s12, s22 }; // Call the step function on the StressUpdateStep class
            // Call the step function on the StressUpdateStep class
            stressStep.stressUpdateHighOrder(
                params, *smesh, stress, { e11, e12, e22 }, hice, cice, deltaT);

            stressDivergence(); // Compute divergence of stress tensor

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

    ModelArray getDG0Data(const std::string& name) override
    {
        if (name == damageName) {
            ModelArray data(ModelArray::Type::H);
            return DGModelArray::dg2ma(damage, data);
        } else {
            return CGDynamicsKernel<DGadvection>::getDG0Data(name);
        }
    }

protected:
    CGVector<CGdegree> avgU;
    CGVector<CGdegree> avgV;

    StressUpdateStep<DGadvection, DGstressDegree>& stressStep;
    const MEBParameters& params;

    Nextsim::DGTransport<DGstressDegree>* stresstransport;

    DGVector<DGadvection> damage;

    // Common brittle parts of the momentum solver.
    void updateMomentum(const TimestepTime& tst) override
    {
#pragma omp parallel for
        for (size_t i = 0; i < u.rows(); ++i) {
            // FIXME dte_over_mass should include snow in the total mass
            const double dteOverMass = deltaT / (params.rho_ice * cgH(i));
            // Memorized initial velocity values
            const double uIce = u(i);
            const double vIce = v(i);

            const double cPrime
                = cgA(i) * params.F_ocean * std::hypot(uOcean(i) - uIce, vOcean(i) - vIce);

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
            const double tauX = dragAtm * uAtmos(i)
                + cPrime * (uOcean(i) * cosOceanAngle - vOcean(i) * sinOceanAngle);
            const double tauY = dragAtm * vAtmos(i)
                + cPrime * (vOcean(i) * cosOceanAngle + uOcean(i) * sinOceanAngle);

            // Stress gradient
            const double gradX = dStressX(i) / pmap->lumpedcgmass(i);
            const double gradY = dStressY(i) / pmap->lumpedcgmass(i);

            u(i) = alpha * uIce + beta * vIce
                + dteOverMass * (alpha * (gradX + tauX) + beta * (gradY + tauY));
            u(i) *= rDenom;

            v(i) = alpha * vIce - beta * uIce
                + dteOverMass * (alpha * (gradY + tauY) + beta * (gradX + tauX));
            v(i) *= rDenom;

            // Calculate the contribution to the average velocity
            avgU(i) += u(i) / nSteps;
            avgV(i) += v(i) / nSteps;
        }
    }
};
}

#endif /* BRITTLECGDYNAMICSKERNEL_HPP */
