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
template <int CG, typename Vec>
static KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, CGDOFS(CG), 1> cgToLocal(
    const Vec& vGlobal, DeviceIndex cgi, DeviceIndex cgShift)
{
    if constexpr (CG == 1) {
        Eigen::Matrix<FloatType, CGDOFS(1), 1> vLocal;
        vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + cgShift),
            vGlobal(cgi + 1 + cgShift);
        return vLocal;
    } else {
        Eigen::Matrix<FloatType, CGDOFS(2), 1> vLocal;
        vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + 2), vGlobal(cgi + cgShift),
            vGlobal(cgi + 1 + cgShift), vGlobal(cgi + 2 + cgShift), vGlobal(cgi + 2 * cgShift),
            vGlobal(cgi + 1 + 2 * cgShift), vGlobal(cgi + 2 + 2 * cgShift);
        return vLocal;
    }
}

/*************************************************************/
template <int DGadvection>
void KokkosMEVPDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    buffers.u0DeviceMut = makeKokkosDeviceView("u0", this->u);
    buffers.v0DeviceMut = makeKokkosDeviceView("v0", this->v);
    buffers.u0Device = buffers.u0DeviceMut;
    buffers.v0Device = buffers.v0DeviceMut;

    std::tie(buffers.s11Host, buffers.s11Device) = makeKokkosDualView("s11", this->s11);
    std::tie(buffers.s12Host, buffers.s12Device) = makeKokkosDualView("s12", this->s12);
    std::tie(buffers.s22Host, buffers.s22Device) = makeKokkosDualView("s22", this->s22);
    std::tie(buffers.e11Host, buffers.e11Device) = makeKokkosDualView("e11", this->e11);
    std::tie(buffers.e12Host, buffers.e12Device) = makeKokkosDualView("e12", this->e12);
    std::tie(buffers.e22Host, buffers.e22Device) = makeKokkosDualView("e22", this->e22);

    std::tie(buffers.hiceHost, buffers.hiceDevice) = makeKokkosDualView("hice", this->hice);
    std::tie(buffers.ciceHost, buffers.ciceDevice) = makeKokkosDualView("cice", this->cice);
    
    // does not depend on the data but it can not be allocated in the constructor
    buffers.PSIAdvectDevice
        = makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection, NGP>, true);
    buffers.PSIStressDevice
        = makeKokkosDeviceView("PSI<DGstress, NGP>", PSI<DGstressComp, NGP>, true);

    // parametric map related
    buffers.lumpedcgmassDevice
        = makeKokkosDeviceView("lumpedcgmass", this->pmap->lumpedcgmass, true);
    buffers.divS1Device = makeKokkosDeviceViewMap("divS1", this->pmap->divS1, true);
    buffers.divS2Device = makeKokkosDeviceViewMap("divS2", this->pmap->divS2, true);
    buffers.divMDevice = makeKokkosDeviceViewMap("divM", this->pmap->divM, true);
    buffers.iMgradXDevice = makeKokkosDeviceViewMap("iMgradX", this->pmap->iMgradX, true);
    buffers.iMgradYDevice = makeKokkosDeviceViewMap("iMgradY", this->pmap->iMgradY, true);
    buffers.iMMDevice = makeKokkosDeviceViewMap("iMM", this->pmap->iMM, true);
    buffers.iMJwPSIDevice = makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, true);

    meshData = std::make_unique_ptr<KokkosMeshData>(*this->smesh);

    // These buffers are only used internally. Thus, synchronisation with CPU only needs to happen
    // to save/load the state. todo: read back buffers if needed in outputs
    Kokkos::deep_copy(buffers.s11Device, buffers.s11Host);
    Kokkos::deep_copy(buffers.s12Device, buffers.s12Host);
    Kokkos::deep_copy(buffers.s22Device, buffers.s22Host);
    Kokkos::deep_copy(buffers.e11Device, buffers.e11Host);
    Kokkos::deep_copy(buffers.e12Device, buffers.e12Host);
    Kokkos::deep_copy(buffers.e22Device, buffers.e22Host);
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
    Kokkos::deep_copy(execSpace, buffers.uDevice, buffers.uHost);
    Kokkos::deep_copy(execSpace, buffers.vDevice, buffers.vHost);
    Kokkos::deep_copy(execSpace, buffers.u0DeviceMut, buffers.uDevice);
    Kokkos::deep_copy(execSpace, buffers.v0DeviceMut, buffers.vDevice);

    Kokkos::deep_copy(execSpace, buffers.uOceanDevice, buffers.uOceanHost);
    Kokkos::deep_copy(execSpace, buffers.vOceanDevice, buffers.vOceanHost);

    Kokkos::deep_copy(execSpace, buffers.uAtmosDevice, buffers.uAtmosHost);
    Kokkos::deep_copy(execSpace, buffers.vAtmosDevice, buffers.vAtmosHost);

    Kokkos::deep_copy(execSpace, buffers.hiceDevice, buffers.hiceHost);
    Kokkos::deep_copy(execSpace, buffers.ciceDevice, buffers.ciceHost);
    Kokkos::deep_copy(execSpace, buffers.cgHDevice, buffers.cgHHost);
    Kokkos::deep_copy(execSpace, buffers.cgADevice, buffers.cgAHost);
    execSpace.fence();
    timerUpload.stop();

    for (size_t mevpstep = 0; mevpstep < this->nSteps; ++mevpstep) {
        timerProj.start();
        projVelocityToStrain(
            buffers, this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
        timerProj.stop();

        timerStress.start();
        updateStressHighOrder(buffers, params, alpha);
        timerStress.stop();

        timerDivergence.start();
        computeStressDivergence(
            buffers, this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
        timerDivergence.stop();

        timerMomentum.start();
        updateMomentumDevice(tst, buffers, params, beta);
        timerMomentum.stop();

        timerBoundary.start();
        applyBoundariesDevice(buffers, this->smesh->nx, this->smesh->ny);
        timerBoundary.stop();
    }

    timerDownload.start();
    Kokkos::deep_copy(execSpace, buffers.uHost, buffers.uDevice);
    Kokkos::deep_copy(execSpace, buffers.vHost, buffers.vDevice);
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
void KokkosMEVPDynamicsKernel<DGadvection>::updateStressHighOrder(
    const KokkosBuffers& _buffers, const VPParameters& _params, FloatType _alpha)
{
    const DeviceIndex n = _buffers.s11Device.extent(0);
    Kokkos::parallel_for(
        "updateStressHighOrder", n, KOKKOS_LAMBDA(const DeviceIndex i) {
            auto s11 = makeEigenMap(_buffers.s11Device);
            auto s12 = makeEigenMap(_buffers.s12Device);
            auto s22 = makeEigenMap(_buffers.s22Device);
            auto e11 = makeEigenMap(_buffers.e11Device);
            auto e12 = makeEigenMap(_buffers.e12Device);
            auto e22 = makeEigenMap(_buffers.e22Device);
            auto hice = makeEigenMap(_buffers.hiceDevice);
            auto cice = makeEigenMap(_buffers.ciceDevice);

            auto PSIAdvect = makeEigenMap(_buffers.PSIAdvectDevice);
            auto PSIStress = makeEigenMap(_buffers.PSIStressDevice);

            auto h_gauss = (hice.row(i) * PSIAdvect).array().max(0.0).matrix();
            auto a_gauss = (cice.row(i) * PSIAdvect).array().max(0.0).min(1.0).matrix();

            EdgeVec P = (_params.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp())
                            .matrix();
            const EdgeVec e11_gauss = e11.row(i) * PSIStress;
            const EdgeVec e12_gauss = e12.row(i) * PSIStress;
            const EdgeVec e22_gauss = e22.row(i) * PSIStress;

            const auto DELTA = (_params.DeltaMin * _params.DeltaMin
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array() + e12_gauss.array().square())
                                   .sqrt()
                                   .matrix();

            const auto map = _buffers.iMJwPSIDevice[i];

            const FloatType alphaInv = 1.0 / _alpha;
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
void KokkosMEVPDynamicsKernel<DGadvection>::updateMomentumDevice(const TimestepTime& tst,
    const KokkosBuffers& _buffers, const VPParameters& _params, FloatType beta)
{
    // Update the velocity
    const FloatType SC = 1.0; ///(1.0-pow(1.0+1.0/beta,-1.0*nSteps));
    const FloatType deltaT = tst.step.seconds();

    //      update by a loop.. implicit parts and h-dependent
    Kokkos::parallel_for(
        "updateMomentum", _buffers.uDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto uOcnRel = _buffers.uOceanDevice(i)
                - _buffers.uDevice(i); // note the reversed sign compared to the v component
            auto vOcnRel = _buffers.vDevice(i) - _buffers.vOceanDevice(i);
            const FloatType absatm
                = Kokkos::sqrt(SQR(_buffers.uAtmosDevice(i)) + SQR(_buffers.vAtmosDevice(i)));
            const FloatType absocn = Kokkos::sqrt(
                SQR(uOcnRel) + SQR(vOcnRel)); // note that the sign of uOcnRel is irrelevant here

            _buffers.uDevice(i) = (1.0
                / (_params.rho_ice * _buffers.cgHDevice(i) / deltaT * (1.0 + beta) // implicit parts
                    + _buffers.cgADevice(i) * _params.F_ocean * absocn) // implicit parts
                * (_params.rho_ice * _buffers.cgHDevice(i) / deltaT
                        * (beta * _buffers.uDevice(i)
                            + _buffers.u0Device(i)) // pseudo - timestepping
                    + _buffers.cgADevice(i)
                        * (_params.F_atm * absatm * _buffers.uAtmosDevice(i) + // atm forcing
                            _params.F_ocean * absocn * SC
                                * _buffers.uOceanDevice(i)) // ocean forcing
                    + _params.rho_ice * _buffers.cgHDevice(i) * _params.fc
                        * vOcnRel // cor + surface
                    + _buffers.dStressXDevice(i) / _buffers.lumpedcgmassDevice(i)));
            _buffers.vDevice(i) = (1.0
                / (_params.rho_ice * _buffers.cgHDevice(i) / deltaT * (1.0 + beta) // implicit parts
                    + _buffers.cgADevice(i) * _params.F_ocean * absocn) // implicit parts
                * (_params.rho_ice * _buffers.cgHDevice(i) / deltaT
                        * (beta * _buffers.vDevice(i)
                            + _buffers.v0Device(i)) // pseudo - timestepping
                    + _buffers.cgADevice(i)
                        * (_params.F_atm * absatm * _buffers.vAtmosDevice(i) + // atm forcing
                            _params.F_ocean * absocn * SC
                                * _buffers.vOceanDevice(i)) // ocean forcing
                    + _params.rho_ice * _buffers.cgHDevice(i) * _params.fc
                        * uOcnRel // here the reversed sign of uOcnRel is used
                    + _buffers.dStressYDevice(i) / _buffers.lumpedcgmassDevice(i)));
        });
}

template class KokkosMEVPDynamicsKernel<1>;
template class KokkosMEVPDynamicsKernel<3>;
template class KokkosMEVPDynamicsKernel<6>;

} // namespace Nextsim