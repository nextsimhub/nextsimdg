#include "include/KokkosBrittleCGDynamicsKernel.hpp"
#include "include/KokkosTimer.hpp"

namespace Nextsim {

template <int DGadvection>
KokkosBrittleCGDynamicsKernel<DGadvection>::KokkosBrittleCGDynamicsKernel(
    const MEBParameters& paramsIn)
    : params(paramsIn)
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

    // Degrees to radians as a hex float
    constexpr FloatType radians = 0x1.1df46a2529d39p-6;
    cosOceanAngle = std::cos(radians * params.ocean_turning_angle);
    sinOceanAngle = std::sin(radians * params.ocean_turning_angle);

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
    //    static KokkosTimer<DETAILED_MEASUREMENTS> timerProj("projGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerStress("stressGPU");
    //    static KokkosTimer<DETAILED_MEASUREMENTS> timerDivergence("divGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerMomentum("momentumGPU");
    //    static KokkosTimer<DETAILED_MEASUREMENTS> timerBoundary("bcGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerUpload("uploadGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerDownload("downloadGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerAdvectStress("advectStressGPU");

    // Let DynamicsKernel handle the advection step
    this->advectionAndLimits(tst);

    // Transport and limits for damage
    const double dt = tst.step.seconds();
    this->dgtransport->step(dt, damage);
    Nextsim::LimitMax(damage, 1.0);
    Nextsim::LimitMin(damage, 1e-12);

    //! Perform transport step for stress
    timerAdvectStress.start();
    stressTransportDevice->prepareAdvection(avgUDevice, avgVDevice);
    stressTransportDevice->step(dt, this->s11Device);
    stressTransportDevice->step(dt, this->s12Device);
    stressTransportDevice->step(dt, this->s22Device);
    timerAdvectStress.stop();

    this->prepareIteration({ { hiceName, this->hice }, { ciceName, this->cice } });

    // The timestep for the brittle solver is the solver subtimestep
    this->deltaT = tst.step.seconds() / this->nSteps;

    timerBBM.start();
    timerUpload.start();
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
    timerUpload.stop();

    for (size_t subStep = 0; subStep < this->nSteps; ++subStep) {
        Base::projectVelocityToStrainDevice(this->uDevice, this->vDevice, this->e11Device,
            this->e12Device, this->e22Device, this->meshData->landMaskDevice, this->iMgradXDevice,
            this->iMgradYDevice, this->iMMDevice, this->smesh->nx, this->smesh->ny,
            this->smesh->CoordinateSystem);

        timerStress.start();
        updateStressHighOrderDevice(this->s11Device, this->s12Device, this->s22Device,
            this->e11Device, this->e12Device, this->e22Device, this->hiceDevice, this->ciceDevice,
            this->damageDevice, this->deltaT);
        timerStress.stop();

        Base::computeStressDivergenceDevice(this->dStressXDevice, this->dStressYDevice,
            this->s11Device, this->s12Device, this->s22Device, this->meshData->landMaskDevice,
            this->divS1Device, this->divS2Device, this->divMDevice, this->meshData->dirichletDevice,
            this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);

        timerMomentum.start();
        updateMomentumDevice(this->uDevice, this->vDevice, this->avgUDevice, this->avgVDevice,
            this->cgHDevice, this->cgADevice, this->uAtmosDevice, this->vAtmosDevice,
            this->uOceanDevice, this->vOceanDevice, this->dStressXDevice, this->dStressYDevice,
            this->lumpedCGMassDevice, this->deltaT, this->params, cosOceanAngle, sinOceanAngle,
            this->nSteps);
        timerMomentum.stop();

        Base::applyBoundariesDevice(this->uDevice, this->vDevice, this->meshData->dirichletDevice,
            this->smesh->nx, this->smesh->ny);
    }

    timerDownload.start();
    Kokkos::deep_copy(execSpace, this->uHost, this->uDevice);
    Kokkos::deep_copy(execSpace, this->vHost, this->vDevice);
    Kokkos::deep_copy(execSpace, this->avgUHost, this->avgUDevice);
    Kokkos::deep_copy(execSpace, this->avgVHost, this->avgVDevice);

    Kokkos::deep_copy(execSpace, this->damageHost, this->damageDevice);

    /*    Kokkos::deep_copy(execSpace, this->s11Host, this->s11Device);
        Kokkos::deep_copy(execSpace, this->s12Host, this->s12Device);
        Kokkos::deep_copy(execSpace, this->s22Host, this->s22Device);*/
    timerDownload.stop();
    timerBBM.stop();

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
        return CGDynamicsKernel<DGadvection>::getDG0Data(name);
    }
}

template <int DGadvection>
ModelArray KokkosBrittleCGDynamicsKernel<DGadvection>::getDGData(const std::string& name) const
{
    if (name == damageName) {
        ModelArray data(ModelArray::Type::DG);
        return DGModelArray::dg2ma(damage, data);
    } else {
        return CGDynamicsKernel<DGadvection>::getDGData(name);
    }
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

template class KokkosBrittleCGDynamicsKernel<1>;
template class KokkosBrittleCGDynamicsKernel<3>;
template class KokkosBrittleCGDynamicsKernel<6>;
}