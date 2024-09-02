#include "include/KokkosBrittleCGDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    KokkosCGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    //! Initialize stress transport
    stressTransport = std::make_unique<Nextsim::DGTransport<DGstressComp>>(*this->smesh);
    stressTransport->settimesteppingscheme("rk2");

    damage.resize_by_mesh(*this->smesh);
    avgU.resize_by_mesh(*this->smesh);
    avgV.resize_by_mesh(*this->smesh);

    // Degrees to radians as a hex float
    constexpr FloatType radians = 0x1.1df46a2529d39p-6;
    cosOceanAngle = std::cos(radians * params.ocean_turning_angle);
    sinOceanAngle = std::sin(radians * params.ocean_turning_angle);

    std::tie(avgUHost, avgUDevice) = makeKokkosDualView("avgU", this->avgU);
    std::tie(avgVHost, avgVDevice) = makeKokkosDualView("avgV", this->avgV);

    std::tie(damageHost, damageDevice) = makeKokkosDualView("damage", this->damage);

    // BBM stress update related
    iMJwPSIAdvectDevice = makeKokkosDeviceViewMap("iMJwPSIAdvect", this->pmap->iMJwPSIAdvect, true);
    std::vector<FloatType> cellSize(this->smesh->nelements);
    for (size_t i = 0; i < this->smesh->nelements; ++i) {
        cellSize[i] = this->smesh->h(i);
    }
    cellSizeDevice = makeKokkosDeviceViewMap("cellSize", cellSize, true);
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::update(const TimestepTime& tst)
{
    // Let DynamicsKernel handle the advection step
    this->advectionAndLimits(tst);

    //! Perform transport step for stress
    stressTransport->prepareAdvection(avgU, avgV);
    stressTransport->step(tst.step.seconds(), this->s11);
    stressTransport->step(tst.step.seconds(), this->s12);
    stressTransport->step(tst.step.seconds(), this->s22);

    // Transport and limits for damage
    this->dgtransport->step(tst.step.seconds(), damage);
    Nextsim::LimitMax(damage, 1.0);
    Nextsim::LimitMin(damage, 1e-12);

    this->prepareIteration({ { hiceName, this->hice }, { ciceName, this->cice } });

    // The timestep for the brittle solver is the solver subtimestep
    this->deltaT = tst.step.seconds() / this->nSteps;

    // explicit execution space enables asynchronous execution
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, avgUDevice, 0.0);
    Kokkos::deep_copy(execSpace, avgVDevice, 0.0);

    Kokkos::deep_copy(execSpace, this->uDevice, this->uHost);
    Kokkos::deep_copy(execSpace, this->vDevice, this->vHost);

    Kokkos::deep_copy(execSpace, this->uOceanDevice, this->uOceanHost);
    Kokkos::deep_copy(execSpace, this->vOceanDevice, this->vOceanHost);

    Kokkos::deep_copy(execSpace, this->uAtmosDevice, this->uAtmosHost);
    Kokkos::deep_copy(execSpace, this->vAtmosDevice, this->vAtmosHost);

    Kokkos::deep_copy(execSpace, this->hiceDevice, this->hiceHost);
    Kokkos::deep_copy(execSpace, this->ciceDevice, this->ciceHost);
    Kokkos::deep_copy(execSpace, this->cgHDevice, this->cgHHost);
    Kokkos::deep_copy(execSpace, this->cgADevice, this->cgAHost);

    Kokkos::deep_copy(execSpace, this->damageDevice, this->damageHost);

    // stress is advected in the brittle model so we have to sync with the host
    Kokkos::deep_copy(execSpace, this->s11Device, this->s11Host);
    Kokkos::deep_copy(execSpace, this->s12Device, this->s12Host);
    Kokkos::deep_copy(execSpace, this->s22Device, this->s22Host);

    for (size_t subStep = 0; subStep < this->nSteps; ++subStep) {
        Base::projectVelocityToStrainDevice(this->uDevice, this->vDevice, this->e11Device,
            this->e12Device, this->e22Device, this->meshData->landMaskDevice, this->iMgradXDevice,
            this->iMgradYDevice, this->iMMDevice, this->smesh->nx, this->smesh->ny,
            this->smesh->CoordinateSystem);

        updateStressHighOrderDevice(this->s11Device, this->s12Device, this->s22Device,
            this->e11Device, this->e12Device, this->e22Device, this->PSIAdvectDevice,
            this->PSIStressDevice, this->hiceDevice, this->ciceDevice, this->damageDevice,
            this->iMJwPSIDevice, this->iMJwPSIAdvectDevice, this->cellSizeDevice, this->deltaT,
            params);

        Base::computeStressDivergenceDevice(this->dStressXDevice, this->dStressYDevice,
            this->s11Device, this->s12Device, this->s22Device, this->meshData->landMaskDevice,
            this->divS1Device, this->divS2Device, this->divMDevice, this->meshData->dirichletDevice,
            this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);

        updateMomentumDevice(this->uDevice, this->vDevice, this->avgUDevice, this->avgVDevice,
            this->cgHDevice, this->cgADevice, this->uAtmosDevice, this->vAtmosDevice,
            this->uOceanDevice, this->vOceanDevice, this->dStressXDevice, this->dStressYDevice,
            this->lumpedCGMassDevice, this->deltaT, this->params, cosOceanAngle, sinOceanAngle,
            this->nSteps);

        Base::applyBoundariesDevice(this->uDevice, this->vDevice, this->meshData->dirichletDevice,
            this->smesh->nx, this->smesh->ny);
    }

    Kokkos::deep_copy(execSpace, this->uHost, this->uDevice);
    Kokkos::deep_copy(execSpace, this->vHost, this->vDevice);
    Kokkos::deep_copy(execSpace, this->avgUHost, this->avgUDevice);
    Kokkos::deep_copy(execSpace, this->avgVHost, this->avgVDevice);

    Kokkos::deep_copy(execSpace, this->damageHost, this->damageDevice);

    Kokkos::deep_copy(execSpace, this->s11Host, this->s11Device);
    Kokkos::deep_copy(execSpace, this->s12Host, this->s12Device);
    Kokkos::deep_copy(execSpace, this->s22Host, this->s22Device);

    // Finally, do the base class update
    DynamicsKernel<DGadvection, DGstressComp>::update(tst);
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::updateStressHighOrderDevice(
    const DeviceViewStress& s11Device, const DeviceViewStress& s12Device,
    const DeviceViewStress& s22Device, const ConstDeviceViewStress& e11Device,
    const ConstDeviceViewStress& e12Device, const ConstDeviceViewStress& e22Device,
    const PSIAdvectView& PSIAdvectDevice, const PSIStressView& PSIStressDevice,
    const ConstDeviceViewAdvect& hiceDevice, const ConstDeviceViewAdvect& ciceDevice,
    const DeviceViewAdvect& damageDevice,
    const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix>& iMJwPSIDevice,
    const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapAdvectMatrix>&
        iMJwPSIAdvectDevice,
    const KokkosDeviceMapView<FloatType>& cellSizeDevice, const FloatType deltaT,
    const MEBParameters& params)
{
    constexpr int NGP = KokkosCGDynamicsKernel<DGadvection>::NGP;
    using EdgeVec = Eigen::Matrix<FloatType, 1, NGP * NGP>;

    Kokkos::parallel_for(
        "updateStressHighOrder", s11Device.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto s11 = makeEigenMap(s11Device);
            auto s12 = makeEigenMap(s12Device);
            auto s22 = makeEigenMap(s22Device);
            auto e11 = makeEigenMap(e11Device);
            auto e12 = makeEigenMap(e12Device);
            auto e22 = makeEigenMap(e22Device);
            auto hice = makeEigenMap(hiceDevice);
            auto cice = makeEigenMap(ciceDevice);
            auto damage = makeEigenMap(damageDevice);

            const auto PSIAdvect = makeEigenMap(PSIAdvectDevice);
            const auto PSIStress = makeEigenMap(PSIStressDevice);

            const EdgeVec hGauss = (hice.row(i) * PSIAdvect).array().max(0.0).matrix();
            const EdgeVec aGauss = (cice.row(i) * PSIAdvect).array().max(0.0).min(1.0).matrix();
            EdgeVec dGauss = (damage.row(i) * PSIAdvect).array().max(1e-12).min(1.0).matrix();

            const EdgeVec e11Gauss = e11.row(i) * PSIStress;
            const EdgeVec e12Gauss = e12.row(i) * PSIStress;
            const EdgeVec e22Gauss = e22.row(i) * PSIStress;

            EdgeVec s11Gauss = s11.row(i) * PSIStress;
            EdgeVec s12Gauss = s12.row(i) * PSIStress;
            EdgeVec s22Gauss = s22.row(i) * PSIStress;

            //! Current normal stress for the evaluation of tildeP (Eqn. 1)
            EdgeVec sigma_n = 0.5 * (s11Gauss.array() + s22Gauss.array());

            //! exp(-C(1-A))
            const EdgeVec expC = (params.compaction_param * (1.0 - aGauss.array())).exp().array();

            // Eqn. 25
            const EdgeVec powalphaexpC
                = (dGauss.array() * expC.array()).pow(params.exponent_relaxation_sigma - 1);
            const EdgeVec time_viscous = params.undamaged_time_relaxation_sigma * powalphaexpC;

            //! BBM  Computing tildeP according to (Eqn. 7b and Eqn. 8)
            // (Eqn. 8)
            const EdgeVec Pmax
                = params.P0 * hGauss.array().pow(params.exponent_compression_factor) * expC.array();

            // (Eqn. 7b) Prepare tildeP
            // tildeP must be capped at 1 to get an elastic response
            // (Eqn. 7b) Select case based on sigma_n
            const EdgeVec tildeP
                = (sigma_n.array() < 0.0)
                      .select((-Pmax.array() / sigma_n.array()).min(1.0).matrix(), 0.);

            // multiplicator
            const EdgeVec multiplicator
                = time_viscous.array() / (time_viscous.array() + (1. - tildeP.array()) * deltaT);

            //! Eqn. 9
            const EdgeVec elasticity
                = hGauss.array() * params.young * dGauss.array() * expC.array();

            // Eqn. 12: first factor on RHS
            /* Stiffness matrix
             * / (K:e)11 \       1     /  1  nu    0  \ / e11 \
             * | (K:e)22 |  =  ------- | nu   1    0  | | e22 |
             * \ (K:e)12 /    1 - nu^2 \  0   0  1-nu / \ e12 /
             */

            const EdgeVec Dunit_factor
                = deltaT * elasticity.array() / (1. - (params.nu0 * params.nu0));

            s11Gauss.array()
                += Dunit_factor.array() * (e11Gauss.array() + params.nu0 * e22Gauss.array());
            s22Gauss.array()
                += Dunit_factor.array() * (params.nu0 * e11Gauss.array() + e22Gauss.array());
            s12Gauss.array() += Dunit_factor.array() * e12Gauss.array() * (1. - params.nu0);

            // //! Implicit part of RHS (Eqn. 33)
            s11Gauss.array() *= multiplicator.array();
            s22Gauss.array() *= multiplicator.array();
            s12Gauss.array() *= multiplicator.array();

            sigma_n = 0.5 * (s11Gauss.array() + s22Gauss.array());
            const EdgeVec tau = (0.25 * (s11Gauss.array() - s22Gauss.array()).square()
                + s12Gauss.array().square())
                                    .sqrt();

            const FloatType scale_coef = std::sqrt(0.1 / cellSizeDevice(i));

            //! Eqn. 22
            const EdgeVec cohesion = params.C_lab * scale_coef * hGauss.array();
            //! Eqn. 30
            const EdgeVec compr_strength = params.compr_strength * scale_coef * hGauss.array();

            // Mohr-Coulomb failure using Mssrs. Plante & Tremblay's formulation
            // sigma_s + tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
            EdgeVec dcrit
                = (tau.array() + params.tan_phi * sigma_n.array() > 0.)
                      .select(
                          cohesion.array() / (tau.array() + params.tan_phi * sigma_n.array()), 1.);

            // Compressive failure using Mssrs. Plante & Tremblay's formulation
            dcrit = (sigma_n.array() < -compr_strength.array())
                        .select(-compr_strength.array() / sigma_n.array(), dcrit);

            // Only damage when we're outside
            dcrit = dcrit.array().min(1.0);

            // Eqn. 29
            const EdgeVec td = cellSizeDevice(i)
                * std::sqrt(2. * (1. + params.nu0) * params.rho_ice) / elasticity.array().sqrt();

            // Update damage
            dGauss.array() -= dGauss.array() * (1. - dcrit.array()) * deltaT / td.array();

            // Relax stress in Gauss points
            s11Gauss.array() -= s11Gauss.array() * (1. - dcrit.array()) * deltaT / td.array();
            s12Gauss.array() -= s12Gauss.array() * (1. - dcrit.array()) * deltaT / td.array();
            s22Gauss.array() -= s22Gauss.array() * (1. - dcrit.array()) * deltaT / td.array();

            // INTEGRATION OF STRESS AND DAMAGE
            // get the inverse of the mass matrix scaled with the test-functions in the gauss
            // points, with the gauss weights and with J. This is a 8 x 9 matrix
            const auto iMJwPSI = iMJwPSIDevice[i];
            s11.row(i) = iMJwPSI * s11Gauss.matrix().transpose();
            s12.row(i) = iMJwPSI * s12Gauss.matrix().transpose();
            s22.row(i) = iMJwPSI * s22Gauss.matrix().transpose();

            damage.row(i) = iMJwPSIAdvectDevice[i] * dGauss.matrix().transpose();
        });
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::updateMomentumDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const DeviceViewCG& avgUDevice, const DeviceViewCG& avgVDevice,
    const ConstDeviceViewCG& cgHDevice, const ConstDeviceViewCG& cgADevice,
    const ConstDeviceViewCG& uAtmosDevice, const ConstDeviceViewCG& vAtmosDevice,
    const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
    const ConstDeviceViewCG& dStressXDevice, const ConstDeviceViewCG& dStressYDevice,
    const ConstDeviceViewCG& lumpedCGMassDevice, const FloatType deltaT,
    const MEBParameters& params, FloatType cosOceanAngle, FloatType sinOceanAngle,
    DeviceIndex nSteps)
{
    Kokkos::parallel_for(
        "updateMomentum", uDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            // FIXME dte_over_mass should include snow in the total mass
            const FloatType dteOverMass = deltaT / (params.rho_ice * cgHDevice(i));
            // Memoized initial velocity values
            const FloatType uIce = uDevice(i);
            const FloatType vIce = vDevice(i);

            const FloatType cPrime = cgADevice(i) * params.F_ocean
                * Kokkos::hypot(uOceanDevice(i) - uIce, vOceanDevice(i) - vIce);

            // FIXME grounding term tauB = cBu[i] / std::hypot(uIce, vIce) + u0
            const FloatType tauB = 0.;
            const FloatType alpha = 1 + dteOverMass * (cPrime * cosOceanAngle + tauB);
            /* FIXME latitude needed for spherical cases
             * const FloatType beta = deltaT * params.fc +
             * dteOverMass * cPrime * std::copysign(sinOceanAngle, lat[i]);
             */
            const FloatType beta = deltaT * params.fc + dteOverMass * cPrime * sinOceanAngle;
            const FloatType rDenom = 1 / (SQR(alpha) + SQR(beta));

            // Atmospheric drag
            const FloatType dragAtm
                = cgADevice(i) * params.F_atm * Kokkos::hypot(uAtmosDevice(i), vAtmosDevice(i));
            const FloatType tauX = dragAtm * uAtmosDevice(i)
                + cPrime * (uOceanDevice(i) * cosOceanAngle - vOceanDevice(i) * sinOceanAngle);
            const FloatType tauY = dragAtm * vAtmosDevice(i)
                + cPrime * (vOceanDevice(i) * cosOceanAngle + uOceanDevice(i) * sinOceanAngle);

            // Stress gradient
            const FloatType gradX = dStressXDevice(i) / lumpedCGMassDevice(i);
            const FloatType gradY = dStressYDevice(i) / lumpedCGMassDevice(i);

            uDevice(i) = alpha * uIce + beta * vIce
                + dteOverMass * (alpha * (gradX + tauX) + beta * (gradY + tauY));
            uDevice(i) *= rDenom;

            vDevice(i) = alpha * vIce - beta * uIce
                + dteOverMass * (alpha * (gradY + tauY) + beta * (gradX + tauX));
            vDevice(i) *= rDenom;

            // Calculate the contribution to the average velocity
            avgUDevice(i) += uDevice(i) / nSteps;
            avgVDevice(i) += vDevice(i) / nSteps;
        });
}
/*
template class KokkosBrittleCGDynamicsKernel<1>;
template class KokkosBrittleCGDynamicsKernel<3>;
template class KokkosBrittleCGDynamicsKernel<6>;*/
// because ParametricMomentumMap<CGdegree>::iMJwPSIAdvect does not properly depend on DGadvection we
// can only build this version since the switch is implemented in compile-time we dont really need
// the other versions anyway
template class KokkosBrittleCGDynamicsKernel<DGCOMP>;

}