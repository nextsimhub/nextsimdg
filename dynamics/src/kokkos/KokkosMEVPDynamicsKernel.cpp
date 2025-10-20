/*!
 * @file KokkosMEVPDynamicsKernel.cpp
 *
 * @date 02 Jun 2025
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosMEVPDynamicsKernel.hpp"
#include "include/KokkosMesh.hpp"
#include "include/KokkosTimer.hpp"
#include <include/constants.hpp>

namespace Nextsim {

/*************************************************************/
template <int DGadvection>
KokkosMEVPDynamicsKernel<DGadvection>::KokkosMEVPDynamicsKernel(const VPParameters& paramsIn)
    : KokkosCGDynamicsKernel<DGadvection>(paramsIn)
    , params(paramsIn)
{
}

/*************************************************************/
template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    KokkosCGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    u0DeviceMut = makeKokkosDeviceView("u0", this->u);
    v0DeviceMut = makeKokkosDeviceView("v0", this->v);
    u0Device = u0DeviceMut;
    v0Device = v0DeviceMut;

    // These buffers are only used internally. Thus, synchronisation with CPU only needs to happen
    // to save/load the state. todo: read back buffers if needed in outputs
    Kokkos::deep_copy(this->s11Device, this->s11Host);
    Kokkos::deep_copy(this->s12Device, this->s12Host);
    Kokkos::deep_copy(this->s22Device, this->s22Host);
    Kokkos::deep_copy(this->e11Device, this->e11Host);
    Kokkos::deep_copy(this->e12Device, this->e12Host);
    Kokkos::deep_copy(this->e22Device, this->e22Host);
}

/*************************************************************/
template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::update(const TimestepTime& tst)
{
    static KokkosTimer<true> timerMevp("mevpGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerProj("projGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerStress("stressGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerDivergence("divGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerMomentum("momentumGPU");
    static KokkosTimer<DETAILED_MEASUREMENTS> timerBoundary("bcGPU");
    static KokkosTimer<true> timerUpload("uploadGPU");
    static KokkosTimer<true> timerDownload("downloadGPU");
    static KokkosTimer<true> timerAdvection("advectionGPU");
    static KokkosTimer<true> timerPrepIt("prepItGPU");

    timerUpload.start();
    // explicit execution space enables asynchronous execution
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, this->uDevice, this->uHost);
    Kokkos::deep_copy(execSpace, this->vDevice, this->vHost);
    Kokkos::deep_copy(execSpace, this->u0DeviceMut, this->uDevice);
    Kokkos::deep_copy(execSpace, this->v0DeviceMut, this->vDevice);

    Kokkos::deep_copy(execSpace, this->uOceanDevice, this->uOceanHost);
    Kokkos::deep_copy(execSpace, this->vOceanDevice, this->vOceanHost);

    Kokkos::deep_copy(execSpace, this->uAtmosDevice, this->uAtmosHost);
    Kokkos::deep_copy(execSpace, this->vAtmosDevice, this->vAtmosHost);

    Kokkos::deep_copy(execSpace, this->hiceDevice, this->hiceHost);
    Kokkos::deep_copy(execSpace, this->ciceDevice, this->ciceHost);
    timerUpload.stop();

    timerAdvection.start();
    this->advectDynamicsFields(tst.step.seconds());
    timerAdvection.stop();

    timerPrepIt.start();
    Base::prepareIterationDevice(this->cgHDevice, this->cgADevice, this->hiceDevice,
        this->ciceDevice, *this->dG2CGAdvectInterpolator);
    this->updateGradientOfSeaSurfaceHeight();
    timerPrepIt.stop();

    timerMevp.start();
    // The critical timestep for the VP solver is the advection timestep
    this->deltaT = tst.step.seconds();

    for (size_t mevpstep = 0; mevpstep < params.nSteps; ++mevpstep) {
        timerProj.start();
        Base::projectVelocityToStrainDevice(this->uDevice, this->vDevice, this->e11Device,
            this->e12Device, this->e22Device, this->meshData->landMaskDevice, this->iMgradXDevice,
            this->iMgradYDevice, this->iMMDevice, this->smesh->nx, this->smesh->ny,
            this->smesh->CoordinateSystem);
        timerProj.stop();

        timerStress.start();
        updateStressHighOrderDevice(this->s11Device, this->s12Device, this->s22Device,
            this->e11Device, this->e12Device, this->e22Device, this->PSIAdvectDevice,
            this->PSIStressDevice, this->hiceDevice, this->ciceDevice,
            this->meshData->landMaskDevice, this->iMJwPSIDevice, params, alpha);
        timerStress.stop();

        timerDivergence.start();
        Base::computeStressDivergenceDevice(this->dStressXDevice, this->dStressYDevice,
            this->s11Device, this->s12Device, this->s22Device, this->meshData->landMaskDevice,
            this->divS1Device, this->divS2Device, this->divMDevice, this->cgLandMaskDevice,
            this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
        timerDivergence.stop();

        timerMomentum.start();
        updateMomentumDevice(this->uDevice, this->vDevice, this->u0Device, this->v0Device,
            this->cgHDevice, this->cgADevice, this->uAtmosDevice, this->vAtmosDevice,
            this->uOceanDevice, this->vOceanDevice, this->dStressXDevice, this->dStressYDevice,
            this->xGradSeaSurfaceHeightDevice, this->yGradSeaSurfaceHeightDevice,
            this->lumpedCGMassDevice, this->cgLandMaskDevice, tst, params, beta);
        timerMomentum.stop();

        timerBoundary.start();
        Base::applyBoundariesDevice(
            this->uDevice, this->vDevice, this->cgLandMaskDevice, this->smesh->nx, this->smesh->ny);
        timerBoundary.stop();
    }
    timerMevp.stop();

    Base::updateIceOceanStressDevice(this->uIceOceanStressDevice, this->vIceOceanStressDevice,
        this->uDevice, this->vDevice, this->uOceanDevice, this->vOceanDevice, this->params,
        this->cosOceanAngle, this->sinOceanAngle);
    Kokkos::deep_copy(execSpace, this->uIceOceanStressHost, this->uIceOceanStressDevice);
    Kokkos::deep_copy(execSpace, this->vIceOceanStressHost, this->vIceOceanStressDevice);

    timerDownload.start();
    Kokkos::deep_copy(execSpace, this->uHost, this->uDevice);
    Kokkos::deep_copy(execSpace, this->vHost, this->vDevice);

    Kokkos::deep_copy(execSpace, this->hiceHost, this->hiceDevice);
    Kokkos::deep_copy(execSpace, this->ciceHost, this->ciceDevice);
    timerDownload.stop();

    // Finally, do the base class update
    DynamicsKernel<DGadvection, DGstressComp>::update(tst);
}

/*************************************************************/
template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::updateStressHighOrderDevice(
    const DeviceViewStress& s11Device, const DeviceViewStress& s12Device,
    const DeviceViewStress& s22Device, const ConstDeviceViewStress& e11Device,
    const ConstDeviceViewStress& e12Device, const ConstDeviceViewStress& e22Device,
    const PSIAdvectView& PSIAdvectDevice, const PSIStressView& PSIStressDevice,
    const ConstDeviceViewAdvect& hiceDevice, const ConstDeviceViewAdvect& ciceDevice,
    const ConstDeviceBitset& landMaskDevice, const GaussMapDevice& iMJwPSIDevice,
    const VPParameters& params, FloatType alpha)
{
    const DeviceIndex n = s11Device.extent(0);
    Kokkos::parallel_for(
        "updateStressHighOrder", n, KOKKOS_LAMBDA(const DeviceIndex i) {
            if (!landMaskDevice.test(i)) {
                return;
            }

            auto s11 = makeEigenMap(s11Device);
            auto s12 = makeEigenMap(s12Device);
            auto s22 = makeEigenMap(s22Device);
            auto e11 = makeEigenMap(e11Device);
            auto e12 = makeEigenMap(e12Device);
            auto e22 = makeEigenMap(e22Device);
            auto hice = makeEigenMap(hiceDevice);
            auto cice = makeEigenMap(ciceDevice);

            const auto PSIAdvect = makeEigenMap(PSIAdvectDevice);
            const auto PSIStress = makeEigenMap(PSIStressDevice);

            auto hGauss = (hice.row(i) * PSIAdvect).array().max(0.0).matrix();
            auto aGauss = (cice.row(i) * PSIAdvect).array().max(0.0).min(1.0).matrix();

            const EdgeVec P = (params.pStar * hGauss.array()
                * (params.compactionParam * (1.0 - aGauss.array())).exp())
                                  .matrix();

            const EdgeVec e11Gauss = e11.row(i) * PSIStress;
            const EdgeVec e12Gauss = e12.row(i) * PSIStress;
            const EdgeVec e22Gauss = e22.row(i) * PSIStress;

            const auto DELTA = (params.deltaMin * params.deltaMin
                + 1.25 * (e11Gauss.array().square() + e22Gauss.array().square())
                + 1.50 * e11Gauss.array() * e22Gauss.array() + e12Gauss.array().square())
                                   .sqrt()
                                   .matrix();

            const auto map = iMJwPSIDevice[i];

            const FloatType alphaInv = 1.0 / alpha;
            const FloatType fac = 1.0 - alphaInv;
            const EdgeVec PDelta = P.array() / DELTA.array();
            s11.row(i) = fac * s11.row(i)
                + (map
                    * (alphaInv
                        * (PDelta.array()
                                * ((5.0 / 8.0) * e11Gauss.array() + (3.0 / 8.0) * e22Gauss.array())
                            - 0.5 * P.array())
                            .matrix()
                            .transpose()))
                      .transpose();
            s22.row(i) = fac * s22.row(i)
                + (map
                    * (alphaInv
                        * (PDelta.array()
                                * ((5.0 / 8.0) * e22Gauss.array() + (3.0 / 8.0) * e11Gauss.array())
                            - 0.5 * P.array())
                            .matrix()
                            .transpose()))
                      .transpose();
            s12.row(i) = fac * s12.row(i)
                + (map
                    * (alphaInv
                        * (PDelta.array() * (1.0 / 4.0) * e12Gauss.array()).matrix().transpose()))
                      .transpose();
        });
}

/*************************************************************/
template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::updateMomentumDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const ConstDeviceViewCG& u0Device,
    const ConstDeviceViewCG& v0Device, const ConstDeviceViewCG& cgHDevice,
    const ConstDeviceViewCG& cgADevice, const ConstDeviceViewCG& uAtmosDevice,
    const ConstDeviceViewCG& vAtmosDevice, const ConstDeviceViewCG& uOceanDevice,
    const ConstDeviceViewCG& vOceanDevice, const ConstDeviceViewCG& dStressXDevice,
    const ConstDeviceViewCG& dStressYDevice, const ConstDeviceViewCG& xGradSeaSurfaceHeightDevice,
    const ConstDeviceViewCG& yGradSeaSurfaceHeight, const ConstDeviceViewCG& lumpedCGMassDevice,
    const ConstDeviceBitset& cgLandMaskDevice, const TimestepTime& tst, const VPParameters& params,
    FloatType beta)
{
    // Update the velocity
    const FloatType deltaT = tst.step.seconds();
    const FloatType FOcean = params.COcean * params.rhoOcean;
    const FloatType FAtm = params.CAtm * params.rhoAtm;

    //      update by a loop.. implicit parts and h-dependent
    Kokkos::parallel_for(
        "updateMomentum", uDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            if (!cgLandMaskDevice.test(i)) {
                return;
            }

            const FloatType u = uDevice(i);
            const FloatType v = vDevice(i);
            const FloatType uOcnRel = u - uOceanDevice(i);
            const FloatType vOcnRel = v - vOceanDevice(i);
            const FloatType absatm = Kokkos::sqrt(sqr(uAtmosDevice(i)) + sqr(vAtmosDevice(i)));
            // note that the sign of uOcnRel is irrelevant here
            const FloatType absocn = Kokkos::sqrt(sqr(uOcnRel) + sqr(vOcnRel));

            // TODO: Take the sign of lat into account for Coriolis term
            uDevice(i) = (1.0
                / (params.rhoIce * cgHDevice(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgADevice(i) * FOcean * absocn) // implicit parts
                * (params.rhoIce * cgHDevice(i) / deltaT
                        * (beta * u + u0Device(i)) // pseudo - timestepping
                    + cgADevice(i)
                        * (FAtm * absatm * uAtmosDevice(i) + // atm forcing
                            FOcean * absocn * uOceanDevice(i)) // ocean forcing
                    + params.rhoIce * cgHDevice(i) * params.fc * v // Coriolis
                    - params.rhoIce * cgHDevice(i) * PhysicalConstants::g
                        * xGradSeaSurfaceHeightDevice(i) // sea surface
                    + dStressXDevice(i) / lumpedCGMassDevice(i))); // internal stress term
            vDevice(i) = (1.0
                / (params.rhoIce * cgHDevice(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgADevice(i) * FOcean * absocn) // implicit parts
                * (params.rhoIce * cgHDevice(i) / deltaT
                        * (beta * v + v0Device(i)) // pseudo - timestepping
                    + cgADevice(i)
                        * (FAtm * absatm * vAtmosDevice(i) + // atm forcing
                            FOcean * absocn * vOceanDevice(i)) // ocean forcing
                    + params.rhoIce * cgHDevice(i) * params.fc * u // Coriolis
                    - params.rhoIce * cgHDevice(i) * PhysicalConstants::g
                        * yGradSeaSurfaceHeight(i) // sea surface
                    + dStressYDevice(i) / lumpedCGMassDevice(i))); // internal stress term
        });
}

template class KokkosMEVPDynamicsKernel<DGCOMP>;

} // namespace Nextsim