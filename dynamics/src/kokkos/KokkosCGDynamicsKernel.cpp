/*!
 * @file KokkosCGDynamicsKernel.cpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosCGDynamicsKernel.hpp"

namespace Nextsim {

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    // velocity
    std::tie(uHost, uDevice) = makeKokkosDualView("u", this->u);
    std::tie(vHost, vDevice) = makeKokkosDualView("v", this->v);

    std::tie(cgHHost, cgHDevice) = makeKokkosDualView("cgH", this->cgH);
    std::tie(cgAHost, cgADevice) = makeKokkosDualView("cgA", this->cgA);

    std::tie(dStressXHost, dStressXDevice) = makeKokkosDualView("dStressX", this->dStressX);
    std::tie(dStressYHost, dStressYDevice) = makeKokkosDualView("dStressY", this->dStressY);

    std::tie(uOceanHost, uOceanDevice) = makeKokkosDualView("uOcean", this->uOcean);
    std::tie(vOceanHost, vOceanDevice) = makeKokkosDualView("vOcean", this->vOcean);

    std::tie(uAtmosHost, uAtmosDevice) = makeKokkosDualView("uAtmos", this->uAtmos);
    std::tie(vAtmosHost, vAtmosDevice) = makeKokkosDualView("vAtmos", this->vAtmos);

    // stress
    std::tie(s11Host, s11Device) = makeKokkosDualView("s11", this->s11);
    std::tie(s12Host, s12Device) = makeKokkosDualView("s12", this->s12);
    std::tie(s22Host, s22Device) = makeKokkosDualView("s22", this->s22);
    std::tie(e11Host, e11Device) = makeKokkosDualView("e11", this->e11);
    std::tie(e12Host, e12Device) = makeKokkosDualView("e12", this->e12);
    std::tie(e22Host, e22Device) = makeKokkosDualView("e22", this->e22);

    assert(this->pmap);
    divS1Device = makeKokkosDeviceViewMap("divS1", this->pmap->divS1, true);
    divS2Device = makeKokkosDeviceViewMap("divS2", this->pmap->divS2, true);
    divMDevice = makeKokkosDeviceViewMap("divM", this->pmap->divM, true);
    iMgradXDevice = makeKokkosDeviceViewMap("iMgradX", this->pmap->iMgradX, true);
    iMgradYDevice = makeKokkosDeviceViewMap("iMgradY", this->pmap->iMgradY, true);
    iMMDevice = makeKokkosDeviceViewMap("iMM", this->pmap->iMM, true);

    assert(this->smesh);
    meshData = std::make_unique<KokkosMeshData>(this->smesh);
}

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

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::projectVelocityToStrainDevice(
    const ConstDeviceViewCG& uDevice, const ConstDeviceViewCG& vDevice,
    const DeviceViewStress& e11Device, const DeviceViewStress& e12Device,
    const DeviceViewStress& e22Device, const ConstDeviceBitset& landMaskDevice,
    const GradMapDevice& iMgradXDevice, const GradMapDevice& iMgradYDevice,
    const GradMapDevice& iMMDevice, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
{
    const DeviceIndex cgshift = CGdegree * nx + 1; //!< Index shift for each row

    // 1D loop is much faster than 2D loop on CPU
    Kokkos::parallel_for(
        "projectVelocityToStrain", nx * ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            const DeviceIndex col = idx % nx;
            const DeviceIndex row = idx / nx;
            const DeviceIndex dgi = nx * row + col; //!< Index of dg vector
            const DeviceIndex cgi
                = CGdegree * cgshift * row + col * CGdegree; //!< Lower left index of cg vector

            // only on ice
            if (!landMaskDevice.test(dgi)) {
                return;
            }

            const auto u = makeEigenMap(uDevice);
            const auto v = makeEigenMap(vDevice);

            auto e11 = makeEigenMap(e11Device);
            auto e12 = makeEigenMap(e12Device);
            auto e22 = makeEigenMap(e22Device);

            // get the local x/y - velocity coefficients on the element
            const auto vx_local = cgToLocal<CGdegree>(u, cgi, cgshift);
            const auto vy_local = cgToLocal<CGdegree>(v, cgi, cgshift);

            // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
            // by integrating rhs and inverting with dG(stress) mass matrix
            const auto iMgradX = iMgradXDevice[dgi];
            e11.row(dgi) = iMgradX * vx_local;
            const auto iMgradY = iMgradYDevice[dgi];
            e22.row(dgi) = iMgradY * vy_local;
            e12.row(dgi) = 0.5 * (iMgradX * vy_local + iMgradY * vx_local);

            if (coordinates == SPHERICAL) {
                const auto iMM = iMMDevice[dgi];
                e11.row(dgi) -= iMM * vy_local;
                e12.row(dgi) += 0.5 * iMM * vx_local;
            }
        });
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::dirichletZero(const DeviceViewCG& v, DeviceIndex nx,
    DeviceIndex ny, const KokkosMeshData::DirichletData& dirichlet)
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
        "dirichletZeroTop", dirichlet[2].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
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

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::computeStressDivergenceDevice(
    const DeviceViewCG& dStressXDevice, const DeviceViewCG& dStressYDevice,
    const ConstDeviceViewStress& s11Device, const ConstDeviceViewStress& s12Device,
    const ConstDeviceViewStress& s22Device, const ConstDeviceBitset& landMaskDevice,
    const DivMapDevice& divS1Device, const DivMapDevice& divS2Device,
    const DivMapDevice& divMDevice, const KokkosMeshData::DirichletData& dirichletDevice,
    DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
{
    using CGVec = Eigen::Vector<Nextsim::FloatType, CGdof>;
    //    static PerfTimer timerDivZero("divZeroGPU");
    // static PerfTimer timerDivComp("divCompGPU");
    //   static PerfTimer timerDivDirichlet("divDirichletGPU");
    // zero buffers
    //   timerDivZero.start();
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, dStressXDevice, 0.0);
    Kokkos::deep_copy(execSpace, dStressYDevice, 0.0);
    //  execSpace.fence();
    //    timerDivZero.stop();

    //   timerDivComp.start();
    Kokkos::parallel_for(
        "computeStressDivergence", nx * ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            const DeviceIndex cx = idx % nx;
            const DeviceIndex cy = idx / nx;
            const DeviceIndex eid = cx + nx * cy;
            // only on ice!
            if (!landMaskDevice.test(eid)) {
                return;
            }

            const auto s11 = makeEigenMap(s11Device);
            const auto s12 = makeEigenMap(s12Device);
            const auto s22 = makeEigenMap(s22Device);

            const auto divS1 = divS1Device[eid];
            const auto divS2 = divS2Device[eid];
            CGVec tx = divS1 * s11.row(eid).transpose() + divS2 * s12.row(eid).transpose();
            CGVec ty = divS1 * s12.row(eid).transpose() + divS2 * s22.row(eid).transpose();

            if (coordinates == SPHERICAL) {
                const auto divM = divMDevice[eid];
                tx += divM * s12.row(eid).transpose();
                ty -= divM * s11.row(eid).transpose();
            }
            const DeviceIndex cgRow = CGdegree * nx + 1;
            //!< lower left CG-index in element (cx,cy)
            const DeviceIndex cg_i = CGdegree * cgRow * cy + CGdegree * cx;

            // Fill the stress divergence values
            for (DeviceIndex row = 0; row <= CGdegree; ++row) {
                for (DeviceIndex col = 0; col <= CGdegree; ++col) {
                    const DeviceIndex dst_idx = cg_i + col + row * cgRow;
                    const DeviceIndex src_idx = col + (CGdegree + 1) * row;
                    Kokkos::atomic_sub(&dStressXDevice(dst_idx), tx(src_idx));
                    Kokkos::atomic_sub(&dStressYDevice(dst_idx), ty(src_idx));
                }
            }
        });
    //    timerDivComp.stop();
    // set zero on the Dirichlet boundaries
    //   timerDivDirichlet.start();
    dirichletZero(dStressXDevice, nx, ny, dirichletDevice);
    dirichletZero(dStressYDevice, nx, ny, dirichletDevice);
    //  timerDivDirichlet.stop();
    static int count = 0;
    ++count;
    if (count % 64 == 0) {
        //    timerDivZero.print();
        //    timerDivComp.print();
        //    timerDivDirichlet.print();
    }
    // todo: add the contributions on the periodic boundaries
    // VectorManipulations::CGAveragePeriodic(*smesh, tx);
    //   VectorManipulations::CGAveragePeriodic(*smesh, ty);
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::applyBoundariesDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const KokkosMeshData::DirichletData& dirichletDevice,
    DeviceIndex nx, DeviceIndex ny)
{
    dirichletZero(uDevice, nx, ny, dirichletDevice);
    dirichletZero(vDevice, nx, ny, dirichletDevice);

    // TODO Periodic boundary conditions.
}

template class KokkosCGDynamicsKernel<1>;
template class KokkosCGDynamicsKernel<3>;
template class KokkosCGDynamicsKernel<6>;

}