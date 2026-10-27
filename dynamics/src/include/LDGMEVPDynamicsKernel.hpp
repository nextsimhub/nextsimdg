/*!
 * @file LDGMEVPDynamicsKernel.hpp
 *
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef LDGMEVPDYNAMICSKERNEL_HPP
#define LDGMEVPDYNAMICSKERNEL_HPP

#include "DynamicsKernel.hpp"
#include "LDGMomentum.hpp"
#include "SlopeLimiter.hpp"
#include "VPParameters.hpp"
#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include <memory>
#include <stdexcept>

namespace Nextsim {

//! Parameters of the LDG discretization
struct LDGParameters {
    FloatType a = 0.4; //!< alternating flux u_hat = (0.5+a) u1 + (0.5-a) u2
    FloatType b = 1.e9; //!< jump penalty b/h
    FloatType alpha = 1500.; //!< mEVP parameters
    FloatType beta = 1500.;
    FloatType pcgTolerance = 1.e-8; //!< relative residual of the PCG solver for the velocity
    int pcgMaxIterations = 1000;
};

/*!
 * The mEVP dynamics with the LDG discretization of the momentum equation: velocity, strain
 * and stress in the DG space with DG components, transport with the full DG velocity.
 *
 * Velocity u = 0 on edges between ice and land, the outer mesh boundary is free. Only
 * CARTESIAN coordinates, no periodic boundaries, no MPI.
 */
template <int DGadvection, int DG>
class LDGMEVPDynamicsKernel : public DynamicsKernel<DGadvection, DG> {
    using Base = DynamicsKernel<DGadvection, DG>;
    using Base::cice;
    using Base::dgtransport;
    using Base::e11;
    using Base::e12;
    using Base::e22;
    using Base::hice;
    using Base::s11;
    using Base::s12;
    using Base::s22;
    using Base::seaSurfaceHeight;
    using Base::smesh;

public:
    LDGMEVPDynamicsKernel(const VPParameters& paramsIn, const LDGParameters& ldgParamsIn)
        : params(paramsIn)
        , ldgParams(ldgParamsIn)
    {
    }

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
#ifdef USE_MPI
        throw std::runtime_error("LDGMEVPDynamics: MPI is not supported");
#endif
        if (isSpherical)
            throw std::runtime_error("LDGMEVPDynamics: only Cartesian coordinates are supported");

        Base::initialise(coords, isSpherical, mask);

        // u = 0 on edges between ice and land, the outer boundary is free
        smesh->dirichletFromLandmask(false);

        momentum = std::make_unique<LDGMomentum<DG, DGadvection>>(*smesh, ldgParams.a, ldgParams.b);
        momentum->pcgTolerance = ldgParams.pcgTolerance;
        momentum->pcgMaxIterations = ldgParams.pcgMaxIterations;

        uIceOceanStress.resize_by_mesh(*smesh);
        vIceOceanStress.resize_by_mesh(*smesh);
        uIceOceanStress.setZero();
        vIceOceanStress.setZero();

        cosOceanAngle = std::cos(radians(params.oceanTurningAngle));
        sinOceanAngle = std::sin(radians(params.oceanTurningAngle));
    }

    /*!
     * Velocity, wind and ocean velocity are given as dG0 (cell averages) and set in the dG0
     * component, all others are handled by DynamicsKernel.
     */
    void setData(const std::string& name, const ModelArray& data) override
    {
        const std::map<std::string, DGVector<DG>*> targetMap = {
            { uName, &momentum->vx },
            { vName, &momentum->vy },
            { uWindName, &momentum->ax },
            { vWindName, &momentum->ay },
            { uOceanName, &momentum->ox },
            { vOceanName, &momentum->oy },
        };
        if (targetMap.count(name)) {
            DGVector<DG>& target = *targetMap.at(name);
            target.setZero();
            DGModelArray::ma2dg(data, target);
        } else {
            Base::setData(name, data);
        }
    }

    //! Velocity and ice-ocean stress as dG0 (cell averages), all others from DynamicsKernel
    ModelArray getDG0Data(const std::string& name) const override
    {
        if (name == uName) {
            ModelArray data(ModelArray::Type::U);
            return DGModelArray::dg2ma(momentum->vx, data);
        } else if (name == vName) {
            ModelArray data(ModelArray::Type::V);
            return DGModelArray::dg2ma(momentum->vy, data);
        } else if (name == uIOStressName) {
            ModelArray data(ModelArray::Type::U);
            return DGModelArray::dg2ma(uIceOceanStress, data);
        } else if (name == vIOStressName) {
            ModelArray data(ModelArray::Type::V);
            return DGModelArray::dg2ma(vIceOceanStress, data);
        } else {
            return Base::getDG0Data(name);
        }
    }

    void prepareAdvection() override
    {
        dgtransport->prepareAdvectionDG(momentum->vx, momentum->vy);
    }

    DGVector<DGadvection>& advectDGVField(FloatType timestep, DGVector<DGadvection>& field,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity()) override
    {
        dgtransport->step(timestep, field);

        //! Slope Limiting
        SlopeLimiter<DGadvection> SL(*smesh);
        bool limitSlope = false;
        if (lowerLimit > -std::numeric_limits<FloatType>::infinity()) {
            SL.limitMin(field, lowerLimit);
            limitSlope = true;
        }
        if (upperLimit < std::numeric_limits<FloatType>::infinity()) {
            SL.limitMax(field, upperLimit);
            limitSlope = true;
        }
        if (limitSlope)
            SL.limit(field);
        return field;
    }

    void update(const TimestepTime& tst) override
    {
        Base::advectDynamicsFields(tst.step.seconds());

        // The critical timestep for the VP solver is the advection timestep
        const FloatType dt = tst.step.seconds();
        momentum->setSeaSurfaceHeight(seaSurfaceHeight);
        momentum->prepareIteration(params, hice, cice, ldgParams.beta, dt);
        for (int mevpstep = 0; mevpstep < params.nSteps; ++mevpstep)
            momentum->mEVPStep(params, ldgParams.alpha, ldgParams.beta, dt);

        // strain and stress for the diagnostics of DynamicsKernel
        e11 = momentum->E11;
        e12 = momentum->E12;
        e22 = momentum->E22;
        s11 = momentum->S11;
        s12 = momentum->S12;
        s22 = momentum->S22;

        updateIceOceanStress();

        Base::update(tst);
    }

protected:
    // The steps of the CG kernels are all part of LDGMomentum::mEVPStep
    void updateMomentum(const TimestepTime& tst) override { }
    void prepareIteration(const typename Base::DataMap& data) override { }
    void projectVelocityToStrain() override { momentum->projectVelocityToStrain(); }
    void stressDivergence() override { }
    void applyBoundaries() override { }

private:
    const VPParameters& params;
    const LDGParameters& ldgParams;

    std::unique_ptr<LDGMomentum<DG, DGadvection>> momentum;

    //! Ice-ocean stress (dG0)
    DGVector<1> uIceOceanStress, vIceOceanStress;
    FloatType cosOceanAngle, sinOceanAngle;

    //! Ice-ocean stress from the cell averages of the velocities, as in CGDynamicsKernel
    void updateIceOceanStress()
    {
        const FloatType FOcean = params.COcean * params.rhoOcean;
#pragma omp parallel for
        for (size_t i = 0; i < smesh->nelements; ++i) {
            if (smesh->landmask[i] == 0) {
                uIceOceanStress(i) = 0;
                vIceOceanStress(i) = 0;
                continue;
            }
            const FloatType uOceanRel = momentum->ox(i, 0) - momentum->vx(i, 0);
            const FloatType vOceanRel = momentum->oy(i, 0) - momentum->vy(i, 0);
            const FloatType cPrime = FOcean * std::hypot(uOceanRel, vOceanRel);
            uIceOceanStress(i) = cPrime * (uOceanRel * cosOceanAngle - vOceanRel * sinOceanAngle);
            vIceOceanStress(i) = cPrime * (vOceanRel * cosOceanAngle + uOceanRel * sinOceanAngle);
        }
    }
};

} /* namespace Nextsim */

#endif /* LDGMEVPDYNAMICSKERNEL_HPP */
