/*!
 * @file KokkosBrittleCGDynamicsKernel.cpp
 *
 * @date 02 Jun 2025
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosBrittleCGDynamicsKernel.hpp"
#include "include/KokkosDGLimit.hpp"
#include "include/KokkosTimer.hpp"
#include <include/constants.hpp>

namespace Nextsim {

template <int DGadvection>
KokkosBrittleCGDynamicsKernel<DGadvection>::KokkosBrittleCGDynamicsKernel(
    const BBMParameters& paramsIn)
    : KokkosCGDynamicsKernel<DGadvection>(paramsIn)
    , params(paramsIn)
{
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    KokkosCGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    //! Initialize stress transport
    cG2DGStressInterpolator
        = std::make_unique<Interpolations::KokkosCG2DGInterpolator<DGstressComp, CGdegree>>(
            *this->smesh);
    stressTransportDevice = std::make_unique<KokkosDGTransport<DGstressComp>>(
        *this->smesh, *this->meshData, *cG2DGStressInterpolator);
    stressTransportDevice->setTimeSteppingScheme(TimeSteppingScheme::RK2);

    damage.resize_by_mesh(*this->smesh);
    avgU.resize_by_mesh(*this->smesh);
    avgV.resize_by_mesh(*this->smesh);
    damage.zero();
    avgU.zero();
    avgV.zero();

    std::tie(avgUHost, avgUDevice) = makeKokkosDualView("avgU", this->avgU);
    std::tie(avgVHost, avgVDevice) = makeKokkosDualView("avgV", this->avgV);

    std::tie(damageHost, damageDevice) = makeKokkosDualView("damage", this->damage);
}

template <typename Mat> void compare(const std::string& name, const Mat& m1, const Mat& m2)
{
    FloatType normRef = m1.norm();
    FloatType normDiff = (m1 - m2).norm();
    std::cout << name << " - abs: " << normDiff << ", rel: " << normDiff / normRef
              << ", norm: " << normRef << std::endl;
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::update(const TimestepTime& tst)
{
    static KokkosTimer<true> timerBBM("bbmGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerProj("projGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerStress("stressGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerDivergence("divGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerMomentum("momentumGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerBoundary("bcGPU");
    static KokkosTimer<true> timerUpload("uploadGPU");
    static KokkosTimer<true> timerDownload("downloadGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerAdvectStress("advectStressGPU");
    static KokkosTimer<true> timerAdvection("advectionGPU");
    static KokkosTimer<true> timerPrepIt("prepItGPU");

    timerUpload.start();
    // explicit execution space enables asynchronous execution
    auto execSpace = Kokkos::DefaultExecutionSpace();

    Kokkos::deep_copy(execSpace, this->uDevice, this->uHost);
    Kokkos::deep_copy(execSpace, this->vDevice, this->vHost);

    Kokkos::deep_copy(execSpace, this->uOceanDevice, this->uOceanHost);
    Kokkos::deep_copy(execSpace, this->vOceanDevice, this->vOceanHost);

    Kokkos::deep_copy(execSpace, this->uAtmosDevice, this->uAtmosHost);
    Kokkos::deep_copy(execSpace, this->vAtmosDevice, this->vAtmosHost);

    Kokkos::deep_copy(execSpace, this->hiceDevice, this->hiceHost);
    Kokkos::deep_copy(execSpace, this->ciceDevice, this->ciceHost);

    Kokkos::deep_copy(execSpace, this->damageDevice, this->damageHost);
    timerUpload.stop();

    const FloatType dt = tst.step.seconds();
    timerAdvection.start();
    this->advectAndLimit(dt);

    // Transport and limits for damage
    this->dGTransportDevice->step(dt, damageDevice);
    limitMax(damageDevice, 1.0);
    limitMin(damageDevice, 1e-12);

    //! Perform transport step for stress
    timerAdvectStress.start();
    stressTransportDevice->prepareAdvection(avgUDevice, avgVDevice);
    stressTransportDevice->step(dt, this->s11Device);
    stressTransportDevice->step(dt, this->s12Device);
    stressTransportDevice->step(dt, this->s22Device);
    timerAdvectStress.stop();
    timerAdvection.stop();

    timerPrepIt.start();
    Base::prepareIterationDevice(this->cgHDevice, this->cgADevice, this->hiceDevice,
        this->ciceDevice, *this->dG2CGAdvectInterpolator);
    this->updateGradientOfSeaSurfaceHeight();
    timerPrepIt.stop();

    timerBBM.start();
    // The timestep for the brittle solver is the solver subtimestep
    this->deltaT = dt / params.nSteps;

    Kokkos::deep_copy(execSpace, avgUDevice, 0.0);
    Kokkos::deep_copy(execSpace, avgVDevice, 0.0);

    for (size_t subStep = 0; subStep < params.nSteps; ++subStep) {
        timerProj.start();
        Base::projectVelocityToStrainDevice(this->uDevice, this->vDevice, this->e11Device,
            this->e12Device, this->e22Device, this->meshData->landMaskDevice, this->iMgradXDevice,
            this->iMgradYDevice, this->iMMDevice, this->smesh->nx, this->smesh->ny,
            this->smesh->CoordinateSystem);
        timerProj.stop();

        timerStress.start();
        updateStressHighOrderDevice(this->s11Device, this->s12Device, this->s22Device,
            this->e11Device, this->e12Device, this->e22Device, this->hiceDevice, this->ciceDevice,
            this->damageDevice, this->deltaT);
        timerStress.stop();

        timerDivergence.start();
        Base::computeStressDivergenceDevice(this->dStressXDevice, this->dStressYDevice,
            this->s11Device, this->s12Device, this->s22Device, this->meshData->landMaskDevice,
            this->divS1Device, this->divS2Device, this->divMDevice, this->cgLandMaskDevice,
            this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
        timerDivergence.stop();

        timerMomentum.start();
        updateMomentumDevice(this->uDevice, this->vDevice, this->avgUDevice, this->avgVDevice,
            this->cgHDevice, this->cgADevice, this->uAtmosDevice, this->vAtmosDevice,
            this->uOceanDevice, this->vOceanDevice, this->dStressXDevice, this->dStressYDevice,
            this->xGradSeaSurfaceHeightDevice, this->yGradSeaSurfaceHeightDevice,
            this->lumpedCGMassDevice, this->cgLandMaskDevice, this->deltaT, this->params,
            this->cosOceanAngle, this->sinOceanAngle, params.nSteps);
        timerMomentum.stop();

        timerBoundary.start();
        Base::applyBoundariesDevice(
            this->uDevice, this->vDevice, this->cgLandMaskDevice, this->smesh->nx, this->smesh->ny);
        timerBoundary.stop();
    }
    timerBBM.stop();

    Base::updateIceOceanStressDevice(this->uIceOceanStressDevice, this->vIceOceanStressDevice,
        this->avgUDevice, this->avgVDevice, this->uOceanDevice, this->vOceanDevice, this->params,
        this->cosOceanAngle, this->sinOceanAngle);
    Kokkos::deep_copy(execSpace, this->uIceOceanStressHost, this->uIceOceanStressDevice);
    Kokkos::deep_copy(execSpace, this->vIceOceanStressHost, this->vIceOceanStressDevice);

    timerDownload.start();
    Kokkos::deep_copy(execSpace, this->uHost, this->uDevice);
    Kokkos::deep_copy(execSpace, this->vHost, this->vDevice);

    Kokkos::deep_copy(execSpace, this->hiceHost, this->hiceDevice);
    Kokkos::deep_copy(execSpace, this->ciceHost, this->ciceDevice);
    Kokkos::deep_copy(execSpace, this->damageHost, this->damageDevice);
    /*    Kokkos::deep_copy(execSpace, this->s11Host, this->s11Device);
        Kokkos::deep_copy(execSpace, this->s12Host, this->s12Device);
        Kokkos::deep_copy(execSpace, this->s22Host, this->s22Device);*/
    timerDownload.stop();

    // Finally, do the base class update
    DynamicsKernel<DGadvection, DGstressComp>::update(tst);
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::setData(
    const std::string& name, const ModelArray& data)
{
    if (name == damageName) {
        DGModelArray::ma2dg(data, damage);
    } else {
        CGDynamicsKernel<DGadvection>::setData(name, data);
    }
}

template <int DGadvection>
ModelArray KokkosBrittleCGDynamicsKernel<DGadvection>::getDG0Data(const std::string& name) const
{
    if (name == damageName) {
        ModelArray data(ModelArray::Type::H);
        return DGModelArray::dg2ma(damage, data);
    } else {
        return KokkosCGDynamicsKernel<DGadvection>::getDG0Data(name);
    }
}

template <int DGadvection>
void KokkosBrittleCGDynamicsKernel<DGadvection>::updateMomentumDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const DeviceViewCG& avgUDevice, const DeviceViewCG& avgVDevice,
    const ConstDeviceViewCG& cgHDevice, const ConstDeviceViewCG& cgADevice,
    const ConstDeviceViewCG& uAtmosDevice, const ConstDeviceViewCG& vAtmosDevice,
    const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
    const ConstDeviceViewCG& dStressXDevice, const ConstDeviceViewCG& dStressYDevice,
    const ConstDeviceViewCG& xGradSeaSurfaceHeightDevice,
    const ConstDeviceViewCG& yGradSeaSurfaceHeightDevice,
    const ConstDeviceViewCG& lumpedCGMassDevice, const ConstDeviceBitset& cgLandMaskDevice,
    const FloatType deltaT, const BBMParameters& params, FloatType cosOceanAngle,
    FloatType sinOceanAngle, DeviceIndex nSteps)
{
    const FloatType FOcean = params.COcean * params.rhoOcean;
    const FloatType FAtm = params.CAtm * params.rhoAtm;

    Kokkos::parallel_for(
        "updateMomentum", uDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            if (!cgLandMaskDevice.test(i)) {
                return;
            }

            // FIXME dte_over_mass should include snow in the total mass
            const FloatType dteOverMass = deltaT / (params.rhoIce * cgHDevice(i));
            // Memoized initial velocity values
            const FloatType uIce = uDevice(i);
            const FloatType vIce = vDevice(i);

            const FloatType cPrime = cgADevice(i) * FOcean
                * Kokkos::hypot(uOceanDevice(i) - uIce, vOceanDevice(i) - vIce);

            // FIXME grounding term tauB = cBu[i] / std::hypot(uIce, vIce) + u0
            const FloatType tauB = 0.;
            const FloatType alpha = 1 + dteOverMass * (cPrime * cosOceanAngle + tauB);
            /* FIXME latitude needed for spherical cases
             * const FloatType beta = deltaT * params.fc +
             * dteOverMass * cPrime * std::copysign(sinOceanAngle, lat[i]);
             */
            const FloatType beta = deltaT * params.fc + dteOverMass * cPrime * sinOceanAngle;
            const FloatType rDenom = 1 / (sqr(alpha) + sqr(beta));

            // Atmospheric drag
            const FloatType dragAtm
                = cgADevice(i) * FAtm * Kokkos::hypot(uAtmosDevice(i), vAtmosDevice(i));
            const FloatType tauX = dragAtm * uAtmosDevice(i)
                + cPrime * (uOceanDevice(i) * cosOceanAngle - vOceanDevice(i) * sinOceanAngle);
            const FloatType tauY = dragAtm * vAtmosDevice(i)
                + cPrime * (vOceanDevice(i) * cosOceanAngle + uOceanDevice(i) * sinOceanAngle);

            // Stress gradient
            const FloatType g = params.rhoIce * cgHDevice(i) * PhysicalConstants::g;
            const FloatType gradX
                = dStressXDevice(i) / lumpedCGMassDevice(i) - g * xGradSeaSurfaceHeightDevice(i);
            const FloatType gradY
                = dStressYDevice(i) / lumpedCGMassDevice(i) - g * yGradSeaSurfaceHeightDevice(i);

            const FloatType uIceNew
                = (alpha * uIce + beta * vIce
                      + dteOverMass * (alpha * (gradX + tauX) + beta * (gradY + tauY)))
                * rDenom;
            uDevice(i) = uIceNew;

            const FloatType vIceNew
                = (alpha * vIce - beta * uIce
                      + dteOverMass * (alpha * (gradY + tauY) + beta * (gradX + tauX)))
                * rDenom;
            vDevice(i) = vIceNew;

            // Calculate the contribution to the average velocity
            avgUDevice(i) += uIceNew / nSteps;
            avgVDevice(i) += vIceNew / nSteps;
        });
}

template class KokkosBrittleCGDynamicsKernel<DGCOMP>;
}