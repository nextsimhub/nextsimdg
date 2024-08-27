/*!
 * @file KokkosMEVPDynamicsKernel.cpp
 *
 * @date Mai 31, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosMEVPDynamicsKernel.hpp"

namespace Nextsim {

/*************************************************************/
class PerfTimer {
public:
    PerfTimer(const std::string& _name)
        : m_name(_name)
        , m_count(0)
        , m_total(0.0)
    {
    }

    void start() { m_start = std::chrono::high_resolution_clock::now(); }
    void stop()
    {
        // ensure that gpu tasks are finished
        Kokkos::fence();
        auto end = std::chrono::high_resolution_clock::now();
        ++m_count;
        m_total += std::chrono::duration<double>(end - m_start).count();
    }
    void print()
    {
        std::cout << m_name << " " << m_total << " " << m_total / m_count << " " << m_count
                  << std::endl;
    }

private:
    std::string m_name;
    int m_count;
    double m_total;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

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

    std::tie(hiceHost, hiceDevice) = makeKokkosDualView("hice", this->hice);
    std::tie(ciceHost, ciceDevice) = makeKokkosDualView("cice", this->cice);

    // does not depend on the data but it can not be allocated in the constructor
    PSIAdvectDevice = makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection, NGP>, true);
    PSIStressDevice = makeKokkosDeviceView("PSI<DGstress, NGP>", PSI<DGstressComp, NGP>, true);

    // parametric map related
    lumpedcgmassDevice = makeKokkosDeviceView("lumpedcgmass", this->pmap->lumpedcgmass, true);
    iMJwPSIDevice = makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, true);

    // These buffers are only used internally. Thus, synchronisation with CPU only needs to happen
    // to save/load the state. todo: read back buffers if needed in outputs
    Kokkos::deep_copy(this->s11Device, this->s11Host);
    Kokkos::deep_copy(this->s12Device, this->s12Host);
    Kokkos::deep_copy(this->s22Device, this->s22Host);
    Kokkos::deep_copy(this->e11Device, this->e11Host);
    Kokkos::deep_copy(this->e12Device, this->e12Host);
    Kokkos::deep_copy(this->e22Device, this->e22Host);
}

template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::update(const TimestepTime& tst)
{
    static PerfTimer timerProj("projGPU");
    static PerfTimer timerStress("stressGPU");
    static PerfTimer timerMevp("mevpGPU");
    static PerfTimer timerDivergence("divGPU");
    static PerfTimer timerMomentum("momentumGPU");
    static PerfTimer timerBoundary("bcGPU");
    static PerfTimer timerUpload("uploadGPU");
    static PerfTimer timerDownload("downloadGPU");

    // Let DynamicsKernel handle the advection step
    DynamicsKernel<DGadvection, DGstressComp>::advectionAndLimits(tst);
    this->prepareIteration({ { hiceName, this->hice }, { ciceName, this->cice } });

    // The critical timestep for the VP solver is the advection timestep
    this->deltaT = tst.step.seconds();

    timerMevp.start();
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
    Kokkos::deep_copy(execSpace, this->cgHDevice, this->cgHHost);
    Kokkos::deep_copy(execSpace, this->cgADevice, this->cgAHost);
    execSpace.fence();
    timerUpload.stop();

    for (size_t mevpstep = 0; mevpstep < this->nSteps; ++mevpstep) {
        timerProj.start();
        Base::projectVelocityToStrainDevice(this->uDevice, this->vDevice, this->e11Device,
            this->e12Device, this->e22Device, this->meshData->landMaskDevice, this->iMgradXDevice,
            this->iMgradYDevice, this->iMMDevice, this->smesh->nx, this->smesh->ny,
            this->smesh->CoordinateSystem);
        timerProj.stop();

        timerStress.start();
        updateStressHighOrder(this->s11Device, this->s12Device, this->s22Device, this->e11Device,
            this->e12Device, this->e22Device, this->PSIAdvectDevice, this->PSIStressDevice,
            this->hiceDevice, this->ciceDevice, this->iMJwPSIDevice, params, alpha);
        timerStress.stop();

        timerDivergence.start();
        Base::computeStressDivergenceDevice(this->dStressXDevice, this->dStressYDevice,
            this->s11Device, this->s12Device, this->s22Device, this->meshData->landMaskDevice,
            this->divS1Device, this->divS2Device, this->divMDevice, this->meshData->dirichletDevice,
            this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
        timerDivergence.stop();

        timerMomentum.start();
        updateMomentumDevice(this->uDevice, this->vDevice, this->u0Device, this->v0Device,
            this->cgHDevice, this->cgADevice, this->uAtmosDevice, this->vAtmosDevice,
            this->uOceanDevice, this->vOceanDevice, this->dStressXDevice, this->dStressYDevice,
            this->lumpedcgmassDevice, tst, params, beta);
        timerMomentum.stop();

        timerBoundary.start();
        Base::applyBoundariesDevice(this->uDevice, this->vDevice, this->meshData->dirichletDevice,
            this->smesh->nx, this->smesh->ny);
        timerBoundary.stop();
    }

    timerDownload.start();
    Kokkos::deep_copy(execSpace, this->uHost, this->uDevice);
    Kokkos::deep_copy(execSpace, this->vHost, this->vDevice);
    execSpace.fence();
    timerDownload.stop();

    timerMevp.stop();
    static int macroStep = 0;
    ++macroStep;
    if (macroStep % 16 == 0) {
        timerMevp.print();
        timerProj.print();
        timerStress.print();
        timerDivergence.print();
        timerMomentum.print();
        timerBoundary.print();
        timerUpload.print();
        timerDownload.print();
    }
    // Finally, do the base class update
    DynamicsKernel<DGadvection, DGstressComp>::update(tst);
}

template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::updateStressHighOrder(const DeviceViewStress& s11Device,
    const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
    const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
    const ConstDeviceViewStress& e22Device, const PSIAdvectView& PSIAdvectDevice,
    const PSIStressView& PSIStressDevice, const ConstDeviceViewAdvect& hiceDevice,
    const ConstDeviceViewAdvect& ciceDevice,
    const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix>& iMJwPSIDevice,
    const VPParameters& params, FloatType alpha)
{
    const DeviceIndex n = s11Device.extent(0);
    Kokkos::parallel_for(
        "updateStressHighOrder", n, KOKKOS_LAMBDA(const DeviceIndex i) {
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

            auto h_gauss = (hice.row(i) * PSIAdvect).array().max(0.0).matrix();
            auto a_gauss = (cice.row(i) * PSIAdvect).array().max(0.0).min(1.0).matrix();

            EdgeVec P = (params.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp())
                            .matrix();
            const EdgeVec e11_gauss = e11.row(i) * PSIStress;
            const EdgeVec e12_gauss = e12.row(i) * PSIStress;
            const EdgeVec e22_gauss = e22.row(i) * PSIStress;

            const auto DELTA = (params.DeltaMin * params.DeltaMin
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array() + e12_gauss.array().square())
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
                                * ((5.0 / 8.0) * e11_gauss.array()
                                    + (3.0 / 8.0) * e22_gauss.array())
                            - 0.5 * P.array())
                              .matrix()
                              .transpose()))
                      .transpose();
            s22.row(i) = fac * s22.row(i)
                + (map
                    * (alphaInv
                        * (PDelta.array()
                                * ((5.0 / 8.0) * e22_gauss.array()
                                    + (3.0 / 8.0) * e11_gauss.array())
                            - 0.5 * P.array())
                              .matrix()
                              .transpose()))
                      .transpose();
            s12.row(i) = fac * s12.row(i)
                + (map
                    * (alphaInv
                        * (PDelta.array() * (1.0 / 4.0) * e12_gauss.array()).matrix().transpose()))
                      .transpose();
        });
}

template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::updateMomentumDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const ConstDeviceViewCG& u0Device,
    const ConstDeviceViewCG& v0Device, const ConstDeviceViewCG& cgHDevice,
    const ConstDeviceViewCG& cgADevice, const ConstDeviceViewCG& uAtmosDevice,
    const ConstDeviceViewCG& vAtmosDevice, const ConstDeviceViewCG& uOceanDevice,
    const ConstDeviceViewCG& vOceanDevice, const ConstDeviceViewCG& dStressXDevice,
    const ConstDeviceViewCG& dStressYDevice, const ConstDeviceViewCG& lumpedcgmassDevice,
    const TimestepTime& tst, const VPParameters& params, FloatType beta)
{
    // Update the velocity
    const FloatType SC = 1.0; ///(1.0-pow(1.0+1.0/beta,-1.0*nSteps));
    const FloatType deltaT = tst.step.seconds();

    //      update by a loop.. implicit parts and h-dependent
    Kokkos::parallel_for(
        "updateMomentum", uDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto uOcnRel = uOceanDevice(i)
                - uDevice(i); // note the reversed sign compared to the v component
            auto vOcnRel = vDevice(i) - vOceanDevice(i);
            const FloatType absatm = Kokkos::sqrt(SQR(uAtmosDevice(i)) + SQR(vAtmosDevice(i)));
            const FloatType absocn = Kokkos::sqrt(
                SQR(uOcnRel) + SQR(vOcnRel)); // note that the sign of uOcnRel is irrelevant here

            uDevice(i) = (1.0
                / (params.rho_ice * cgHDevice(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgADevice(i) * params.F_ocean * absocn) // implicit parts
                * (params.rho_ice * cgHDevice(i) / deltaT
                        * (beta * uDevice(i) + u0Device(i)) // pseudo - timestepping
                    + cgADevice(i)
                        * (params.F_atm * absatm * uAtmosDevice(i) + // atm forcing
                            params.F_ocean * absocn * SC * uOceanDevice(i)) // ocean forcing
                    + params.rho_ice * cgHDevice(i) * params.fc * vOcnRel // cor + surface
                    + dStressXDevice(i) / lumpedcgmassDevice(i)));
            vDevice(i) = (1.0
                / (params.rho_ice * cgHDevice(i) / deltaT * (1.0 + beta) // implicit parts
                    + cgADevice(i) * params.F_ocean * absocn) // implicit parts
                * (params.rho_ice * cgHDevice(i) / deltaT
                        * (beta * vDevice(i) + v0Device(i)) // pseudo - timestepping
                    + cgADevice(i)
                        * (params.F_atm * absatm * vAtmosDevice(i) + // atm forcing
                            params.F_ocean * absocn * SC * vOceanDevice(i)) // ocean forcing
                    + params.rho_ice * cgHDevice(i) * params.fc
                        * uOcnRel // here the reversed sign of uOcnRel is used
                    + dStressYDevice(i) / lumpedcgmassDevice(i)));
        });
}

template class KokkosMEVPDynamicsKernel<1>;
template class KokkosMEVPDynamicsKernel<3>;
template class KokkosMEVPDynamicsKernel<6>;

} // namespace Nextsim