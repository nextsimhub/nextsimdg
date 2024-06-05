/*!
 * @file KokkosVPCGDynamicsKernel.hpp
 *
 * @date Mai 31, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "KokkosVPCGDynamicsKernel.hpp"

namespace Nextsim {

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
        auto end = std::chrono::high_resolution_clock::now();
        ++m_count;
        m_total += std::chrono::duration<double>(end - m_start).count();
    }
    void print()
    {
        std::cout << m_name << " total: " << m_total << ", avg: " << m_total / m_count << "\n";
    }

private:
    std::string m_name;
    int m_count;
    double m_total;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

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

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    std::tie(buffers.uHost, buffers.uDevice) = makeKokkosDualView("u", this->u);
    std::tie(buffers.vHost, buffers.vDevice) = makeKokkosDualView("v", this->v);
    buffers.u0Device = makeKokkosDeviceView("u0", this->u);
    buffers.v0Device = makeKokkosDeviceView("v0", this->v);

    std::tie(buffers.dStressXHost, buffers.dStressXDevice)
        = makeKokkosDualView("dStressX", this->dStressX);
    std::tie(buffers.dStressYHost, buffers.dStressYDevice)
        = makeKokkosDualView("dStressY", this->dStressY);

    std::tie(buffers.s11Host, buffers.s11Device) = makeKokkosDualView("s11", this->s11);
    std::tie(buffers.s12Host, buffers.s12Device) = makeKokkosDualView("s12", this->s12);
    std::tie(buffers.s22Host, buffers.s22Device) = makeKokkosDualView("s11", this->s22);
    std::tie(buffers.e11Host, buffers.e11Device) = makeKokkosDualView("e11", this->e11);
    std::tie(buffers.e12Host, buffers.e12Device) = makeKokkosDualView("e12", this->e12);
    std::tie(buffers.e22Host, buffers.e22Device) = makeKokkosDualView("e11", this->e22);

    std::tie(buffers.hiceHost, buffers.hiceDevice) = makeKokkosDualView("hice", this->hice);
    std::tie(buffers.ciceHost, buffers.ciceDevice) = makeKokkosDualView("cice", this->cice);
    std::tie(buffers.cgHHost, buffers.cgHDevice) = makeKokkosDualView("cgH", this->cgH);
    std::tie(buffers.cgAHost, buffers.cgADevice) = makeKokkosDualView("cgA", this->cgA);

    std::tie(buffers.uOceanHost, buffers.uOceanDevice) = makeKokkosDualView("uOcean", this->uOcean);
    std::tie(buffers.vOceanHost, buffers.vOceanDevice) = makeKokkosDualView("vOcean", this->vOcean);

    std::tie(buffers.uAtmosHost, buffers.uAtmosDevice) = makeKokkosDualView("uAtmos", this->uAtmos);
    std::tie(buffers.vAtmosHost, buffers.vAtmosDevice) = makeKokkosDualView("vAtmos", this->vAtmos);

    // does not depend on the data but it can not be allocated in the constructor
    buffers.PSIAdvectDevice
        = makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection, NGP>, true);
    buffers.PSIStressDevice
        = makeKokkosDeviceView("PSI<DGstress, NGP>", PSI<DGstressDegree, NGP>, true);

    // parametric map
    buffers.lumpedcgmassDevice
        = makeKokkosDeviceView("lumpedcgmass", this->pmap->lumpedcgmass, true);
    buffers.divS1Device = makeKokkosDeviceViewMap("divS1", this->pmap->divS1, true);
    buffers.divS2Device = makeKokkosDeviceViewMap("divS2", this->pmap->divS2, true);
    buffers.divMDevice = makeKokkosDeviceViewMap("divM", this->pmap->divM, true);
    buffers.iMgradXDevice = makeKokkosDeviceViewMap("iMgradX", this->pmap->iMgradX, true);
    buffers.iMgradYDevice = makeKokkosDeviceViewMap("iMgradY", this->pmap->iMgradY, true);
    buffers.iMMDevice = makeKokkosDeviceViewMap("iMM", this->pmap->iMM, true);
    buffers.iMJwPSIDevice = makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, true);
    stressStep.setPMap(this->pmap); // only needed for debugging

    // mesh related
    // boundary data
    for (size_t i = 0; i < 4; ++i) {
        buffers.dirichletDevice[i] = makeKokkosDeviceViewMap(
            "dirichlet" + std::to_string(i), this->smesh->dirichlet[i], true);
        //    buffers.periodicDevice[i] = makeKokkosDeviceViewMap(
        //        "periodic" + std::to_string(i), this->smesh->periodic[i], true);
    }
    std::cout << "periodic bcs: " << this->smesh->periodic.size() << "\n";
    for (const auto& bc : this->smesh->periodic) {
        std::cout << ", " << bc.size();
    }
    std::cout << "\n";

    // unfortunately there is no more direct way to initialize an Kokkos::Bitset
    const unsigned nBits = this->smesh->landmask.size();
    Kokkos::Bitset<Kokkos::HostSpace> landMaskHost(nBits);
    landMaskHost.clear();
    for (unsigned i = 0; i < nBits; ++i) {
        if (this->smesh->landmask[i])
            landMaskHost.set(i);
    }
    Kokkos::Bitset<Kokkos::DefaultExecutionSpace> landMaskTemp(nBits);
    Kokkos::deep_copy(landMaskTemp, landMaskHost);
    buffers.landMaskDevice = landMaskTemp;
}

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::update(const TimestepTime& tst)
{
    // Let DynamicsKernel handle the advection step
    DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits(tst);
    this->prepareIteration({ { hiceName, this->hice }, { ciceName, this->cice } });

    // The critical timestep for the VP solver is the advection timestep
    this->deltaT = tst.step.seconds();

    Kokkos::deep_copy(buffers.uDevice, buffers.uHost);
    Kokkos::deep_copy(buffers.vDevice, buffers.vHost);
    Kokkos::deep_copy(buffers.u0Device, buffers.uDevice);
    Kokkos::deep_copy(buffers.v0Device, buffers.vDevice);
    u0 = this->u;
    v0 = this->v;

    Kokkos::deep_copy(buffers.uOceanDevice, buffers.uOceanHost);
    Kokkos::deep_copy(buffers.vOceanDevice, buffers.vOceanHost);

    Kokkos::deep_copy(buffers.uAtmosDevice, buffers.uAtmosHost);
    Kokkos::deep_copy(buffers.vAtmosDevice, buffers.vAtmosHost);

    Kokkos::deep_copy(buffers.s11Device, buffers.s11Host);
    Kokkos::deep_copy(buffers.s12Device, buffers.s12Host);
    Kokkos::deep_copy(buffers.s22Device, buffers.s22Host);
    /*     Kokkos::deep_copy(buffers.e11Device, buffers.e11Host);
         Kokkos::deep_copy(buffers.e12Device, buffers.e12Host);
         Kokkos::deep_copy(buffers.e22Device, buffers.e22Host);*/
    Kokkos::deep_copy(buffers.hiceDevice, buffers.hiceHost);
    Kokkos::deep_copy(buffers.ciceDevice, buffers.ciceHost);
    Kokkos::deep_copy(buffers.cgHDevice, buffers.cgHHost);
    Kokkos::deep_copy(buffers.cgADevice, buffers.cgAHost);

    auto checkDiff = [](const auto& a, const auto& b, bool detailed) {
        if (detailed) {
            std::cout << a.row(7) - b.row(7) << "\n";
        }
        auto aNorm = a.norm();
        std::cout << "a: " << aNorm << ", b: " << b.norm()
                  << ", (a-b)/|a|: " << (a - b).norm() / aNorm << std::endl;
    };
    PerfTimer timerProjGPU("projGPU");
    PerfTimer timerProjCPU("projCPU");

    for (size_t mevpstep = 0; mevpstep < this->nSteps; ++mevpstep) {
        Kokkos::deep_copy(buffers.uDevice, buffers.uHost);
        Kokkos::deep_copy(buffers.vDevice, buffers.vHost);

        Kokkos::deep_copy(buffers.s11Device, buffers.s11Host);
        Kokkos::deep_copy(buffers.s12Device, buffers.s12Host);
        Kokkos::deep_copy(buffers.s22Device, buffers.s22Host);
        Kokkos::deep_copy(buffers.hiceDevice, buffers.hiceHost);
        Kokkos::deep_copy(buffers.ciceDevice, buffers.ciceHost);
        timerProjGPU.start();
        projVelocityToStrain(
            buffers, this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
        Kokkos::fence();
        timerProjGPU.stop();

        timerProjCPU.start();
        this->projectVelocityToStrain();
        timerProjCPU.stop();
        /*      auto tempE11 = this->e11;
              auto tempE12 = this->e12;
              auto tempE22 = this->e22;
              Kokkos::deep_copy(buffers.e11Host, buffers.e11Device);
              Kokkos::deep_copy(buffers.e12Host, buffers.e12Device);
              Kokkos::deep_copy(buffers.e22Host, buffers.e22Device);
        checkDiff(tempE11, this->e11, false);
        checkDiff(tempE12, this->e12, false);
        checkDiff(tempE22, this->e22, false);*/

        /*            std::array<DGVector<DGstressDegree>, N_TENSOR_ELEMENTS> stressTemp {
           this->s11, this->s12, this->s22 };
                    std::array<std::reference_wrapper<DGVector<DGstressDegree>>,
           N_TENSOR_ELEMENTS> stress = { stressTemp[0], stressTemp[1], stressTemp[2] };*/
        stressUpdateHighOrder(buffers, params, alpha);

        Kokkos::deep_copy(buffers.s11Host, buffers.s11Device);
        Kokkos::deep_copy(buffers.s12Host, buffers.s12Device);
        Kokkos::deep_copy(buffers.s22Host, buffers.s22Device);

        this->stressDivergence(); // Compute divergence of stress tensor
        auto tempDStressX = this->dStressX;
        auto tempDStressY = this->dStressY;
        computeStressDivergence(
            buffers, this->smesh->nx, this->smesh->ny, this->smesh->CoordinateSystem);
/*        Kokkos::deep_copy(buffers.dStressXHost, buffers.dStressXDevice);
        Kokkos::deep_copy(buffers.dStressYHost, buffers.dStressYDevice);
        checkDiff(tempDStressX, this->dStressX, false);
        checkDiff(tempDStressY, this->dStressY, false);*/

        updateMomentum(tst);
        auto tempU = this->u;
        auto tempV = this->v;
        updateMomentumDevice(tst, buffers, params, beta);
        Kokkos::deep_copy(buffers.uHost, buffers.uDevice);
        Kokkos::deep_copy(buffers.vHost, buffers.vDevice);
        checkDiff(tempU, this->u, false);
        checkDiff(tempV, this->v, false);
        
        applyBoundariesDevice(buffers, this->smesh->nx, this->smesh->ny);
    }
    //    timerProjGPU.print();
    //    timerProjCPU.print();
    // Finally, do the base class update
    DynamicsKernel<DGadvection, DGstressDegree>::update(tst);
}

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::projVelocityToStrain(
    const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
{
    const DeviceIndex cgshift = CGdegree * nx + 1; //!< Index shift for each row

    // parallelize over 2D grid
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({ 0, 0 }, { nx, ny });
    Kokkos::parallel_for(
        "projectVelocityToStrain", policy,
        KOKKOS_LAMBDA(const DeviceIndex col, const DeviceIndex row) {
            auto u = makeEigenMap(_buffers.uDevice);
            auto v = makeEigenMap(_buffers.vDevice);

            auto e11 = makeEigenMap(_buffers.e11Device);
            auto e12 = makeEigenMap(_buffers.e12Device);
            auto e22 = makeEigenMap(_buffers.e22Device);

            const DeviceIndex dgi = nx * row + col; //!< Index of dg vector
            const DeviceIndex cgi
                = CGdegree * cgshift * row + col * CGdegree; //!< Lower left index of cg vector

            // only on ice
            if (!_buffers.landMaskDevice.test(dgi)) {
                return;
            }

            // get the local x/y - velocity coefficients on the element
            auto vx_local = cgToLocal<CGdegree>(u, cgi, cgshift);
            auto vy_local = cgToLocal<CGdegree>(v, cgi, cgshift);

            // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
            // by integrating rhs and inverting with dG(stress) mass matrix
            e11.row(dgi) = _buffers.iMgradXDevice[dgi] * vx_local;
            //        for(int i = 0; i < 8; ++i)
            //            e11.row(dgi)(i) = _buffers.iMgradXDevice[dgi](i,0);
            e22.row(dgi) = _buffers.iMgradYDevice[dgi] * vy_local;
            e12.row(dgi) = 0.5
                * (_buffers.iMgradXDevice[dgi] * vy_local + _buffers.iMgradYDevice[dgi] * vx_local);

            if (coordinates == SPHERICAL) {
                e11.row(dgi) -= _buffers.iMMDevice[dgi] * vy_local;
                e12.row(dgi) += 0.5 * _buffers.iMMDevice[dgi] * vx_local;
            }
        });
}

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::stressUpdateHighOrder(
    const KokkosBuffers& _buffers, const VPParameters& _params, FloatType _alpha)
{
    const DeviceIndex n = _buffers.s11Device.extent(0);
    Kokkos::parallel_for(
        "stressUpdateHighOrder", n, KOKKOS_LAMBDA(const DeviceIndex i) {
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
/*
template <int DGadvection>
KOKKOS_INLINE_FUNCTION static void addStressTensorCell(const size_t eid, const size_t cx, const
size_t cy)
{
    auto divS1 = _buffers.divS1Device[eid];
    auto divS2 = _buffers.divS2Device[eid];
    Eigen::Vector<Nextsim::FloatType, CGdof> tx
        = (divS1 * s11.row(eid).transpose() + divS2 * s12.row(eid).transpose());
    Eigen::Vector<Nextsim::FloatType, CGdof> ty
        = (divS1 * s12.row(eid).transpose() + divS2 * s22.row(eid).transpose());

    if (coordinates == SPHERICAL) {
        auto divM = _buffers.divMDevice[eid];
        tx += divM * s12.row(eid).transpose();
        ty -= divM * s11.row(eid).transpose();
    }
    const unsigned cgRow = CGdegree * nx + 1;
    const unsigned cg_i
        = CGdegree * cgRow * cy + CGdegree * cx; //!< lower left CG-index in element (cx,cy)

    // Fill the stress divergence values
    for (int row = 0; row <= CGdegree; ++row) {
        for (int col = 0; col <= CGdegree; ++col) {
            _buffers.dStressXDevice(cg_i + col + row * cgRow, 0) -= tx(col + (CGdegree + 1) * row);
            _buffers.dStressYDevice(cg_i + col + row * cgRow, 0) -= ty(col + (CGdegree + 1) * row);
        }
    }
}*/

// todo: should be KokkosVPCGDynamicsKernel::DeviceViewCG
template <typename DeviceViewCG, typename KokkosDeviceMapView>
void dirichletZero(DeviceViewCG& v, DeviceIndex nx, DeviceIndex ny,
    const std::array<KokkosDeviceMapView, 4>& dirichlet)
{
    // bot
    Kokkos::parallel_for(
        "dirichletZeroBot", dirichlet[0].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = dirichlet[0][i];
            const DeviceIndex ix = eid % nx; // compute coordinates of element
            const DeviceIndex iy = eid / nx;
            for (DeviceIndex j = 0; j < CGdegree + 1; ++j) {
                v(iy * CGdegree * (CGdegree * nx + 1) + CGdegree * ix + j) = 0.0;
            }
        });
    // right
    Kokkos::parallel_for(
        "dirichletZeroRight", dirichlet[1].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = dirichlet[1][i];
            const DeviceIndex ix = eid % nx; // compute coordinates of element
            const DeviceIndex iy = eid / nx;
            for (DeviceIndex j = 0; j < CGdegree + 1; ++j) {
                v(iy * CGdegree * (CGdegree * nx + 1) + CGdegree * ix + CGdegree
                    + (CGdegree * nx + 1) * j)
                    = 0.0;
            }
        });
    // top
    Kokkos::parallel_for(
        "dirichletZeroRight", dirichlet[2].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = dirichlet[2][i];
            const DeviceIndex ix = eid % nx; // compute coordinates of element
            const DeviceIndex iy = eid / nx;
            for (DeviceIndex j = 0; j < CGdegree + 1; ++j) {
                v((iy + 1) * CGdegree * (CGdegree * nx + 1) + CGdegree * ix + j) = 0.0;
            }
        });
    // left
    Kokkos::parallel_for(
        "dirichletZeroLeft", dirichlet[3].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = dirichlet[3][i];
            const DeviceIndex ix = eid % nx; // compute coordinates of element
            const DeviceIndex iy = eid / nx;
            for (DeviceIndex j = 0; j < CGdegree + 1; ++j) {
                v(iy * CGdegree * (CGdegree * nx + 1) + CGdegree * ix + (CGdegree * nx + 1) * j)
                    = 0.0;
            }
        });
}
/*
template <typename DeviceViewCG, typename KokkosDeviceMapView>
void CGAveragePeriodic(DeviceViewCG& v, DeviceIndex nx, const std::array<KokkosDeviceMapView, 4>&
periodic)
{
    // the two segments bottom, right, top, left, are each processed in parallel
    for (size_t seg = 0; seg < smesh.periodic.size(); ++seg) {
        // #pragma omp parallel for
        for (size_t i = 0; i < smesh.periodic[seg].size(); ++i) {

            const size_t ptype = smesh.periodic[seg][i][0];
            const size_t eid_lb = smesh.periodic[seg][i][2];
            const size_t eid_rt = smesh.periodic[seg][i][1];

            size_t ix_lb = eid_lb % smesh.nx;
            size_t iy_lb = eid_lb / smesh.nx;
            size_t i0_lb = (CG * smesh.nx + 1) * CG * iy_lb
                + CG * ix_lb; // lower/left index in left/bottom element
            size_t ix_rt = eid_rt % smesh.nx;
            size_t iy_rt = eid_rt / smesh.nx;
            size_t i0_rt = (CG * smesh.nx + 1) * CG * iy_rt
                + CG * ix_rt; // lower/left index in right/top element

            if (ptype == 0) // X-edge, bottom/top
            {
                for (size_t j = 0; j <= CG; ++j) {
                    v(i0_lb + j) = 0.5 * (v(i0_lb + j) + v(i0_rt + CG * (CG * smesh.nx + 1) + j));
                    v(i0_rt + CG * (CG * smesh.nx + 1) + j) = v(i0_lb + j);
                }
            } else if (ptype == 1) // Y-edge, left/right
            {
                for (size_t j = 0; j <= CG; ++j) {
                    const size_t i1 = i0_lb + j * (CG * smesh.nx + 1);
                    const size_t i2 = i0_rt + CG + j * (CG * smesh.nx + 1);
                    v(i1) = 0.5 * (v(i1) + v(i2));
                    v(i2) = v(i1);
                }
            } else
                abort();
        }
    }
    }*/

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::computeStressDivergence(
    const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
{
    // zero buffers
    Kokkos::parallel_for(
        "initStressDivergence", _buffers.dStressXDevice.extent(0),
        KOKKOS_LAMBDA(const DeviceIndex i) {
            _buffers.dStressXDevice(i) = 0.0;
            _buffers.dStressYDevice(i) = 0.0;
        });

    // parallelization in stripes
    for (DeviceIndex p = 0; p < 2; ++p) {
        Kokkos::parallel_for(
            "computeStressDivergence", ny, KOKKOS_LAMBDA(const DeviceIndex cy) {
                auto s11 = makeEigenMap(_buffers.s11Device);
                auto s12 = makeEigenMap(_buffers.s12Device);
                auto s22 = makeEigenMap(_buffers.s22Device);

                //!< loop over all cells of the mesh
                if (cy % 2 == p) {
                    DeviceIndex eid = nx * cy;
                    for (DeviceIndex cx = 0; cx < nx; ++cx, ++eid) {
                        //!< loop over all cells of the mesh
                        // only on ice!
                        if (_buffers.landMaskDevice.test(eid)) {
                            auto divS1 = _buffers.divS1Device[eid];
                            auto divS2 = _buffers.divS2Device[eid];
                            Eigen::Vector<Nextsim::FloatType, CGdof> tx
                                = (divS1 * s11.row(eid).transpose()
                                    + divS2 * s12.row(eid).transpose());
                            Eigen::Vector<Nextsim::FloatType, CGdof> ty
                                = (divS1 * s12.row(eid).transpose()
                                    + divS2 * s22.row(eid).transpose());

                            if (coordinates == SPHERICAL) {
                                auto divM = _buffers.divMDevice[eid];
                                tx += divM * s12.row(eid).transpose();
                                ty -= divM * s11.row(eid).transpose();
                            }
                            const DeviceIndex cgRow = CGdegree * nx + 1;
                            const DeviceIndex cg_i = CGdegree * cgRow * cy
                                + CGdegree * cx; //!< lower left CG-index in element (cx,cy)

                            // Fill the stress divergence values
                            for (DeviceIndex row = 0; row <= CGdegree; ++row) {
                                for (DeviceIndex col = 0; col <= CGdegree; ++col) {
                                    _buffers.dStressXDevice(cg_i + col + row * cgRow)
                                        -= tx(col + (CGdegree + 1) * row);
                                    _buffers.dStressYDevice(cg_i + col + row * cgRow)
                                        -= ty(col + (CGdegree + 1) * row);
                                }
                            }
                        }
                    }
                }
            });
    }
    // set zero on the Dirichlet boundaries
    dirichletZero(_buffers.dStressXDevice, nx, ny, _buffers.dirichletDevice);
    dirichletZero(_buffers.dStressYDevice, nx, ny, _buffers.dirichletDevice);
    // todo: add
    // add the contributions on the periodic boundaries
    // VectorManipulations::CGAveragePeriodic(*smesh, tx);
    //   VectorManipulations::CGAveragePeriodic(*smesh, ty);
}

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::applyBoundariesDevice(
    const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny)
{
    dirichletZero(_buffers.uDevice, nx, ny, _buffers.dirichletDevice);
    dirichletZero(_buffers.vDevice, nx, ny, _buffers.dirichletDevice);
    ;
    // TODO Periodic boundary conditions.
}

template <int DGadvection>
void KokkosVPCGDynamicsKernel<DGadvection>::updateMomentumDevice(const TimestepTime& tst,
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

template class KokkosVPCGDynamicsKernel<1>;
template class KokkosVPCGDynamicsKernel<3>;
template class KokkosVPCGDynamicsKernel<6>;

} // namespace Nextsim