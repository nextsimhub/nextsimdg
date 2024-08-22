/*!
 * @file KokkosCGDynamicsKernel.cpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosCGDynamicsKernel.hpp"

namespace Nextsim {

void KokkosCGDynamicsKernel<DGadvection>::initialise()
{
    CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

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
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::projectVelocityToStrainDevice(
    const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
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
            if (!_buffers.landMaskDevice.test(dgi)) {
                return;
            }

            const auto u = makeEigenMap(_buffers.uDevice);
            const auto v = makeEigenMap(_buffers.vDevice);

            auto e11 = makeEigenMap(_buffers.e11Device);
            auto e12 = makeEigenMap(_buffers.e12Device);
            auto e22 = makeEigenMap(_buffers.e22Device);

            // get the local x/y - velocity coefficients on the element
            const auto vx_local = cgToLocal<CGdegree>(u, cgi, cgshift);
            const auto vy_local = cgToLocal<CGdegree>(v, cgi, cgshift);

            // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
            // by integrating rhs and inverting with dG(stress) mass matrix
            const auto iMgradX = _buffers.iMgradXDevice[dgi];
            e11.row(dgi) = iMgradX * vx_local;
            const auto iMgradY = _buffers.iMgradYDevice[dgi];
            e22.row(dgi) = iMgradY * vy_local;
            e12.row(dgi) = 0.5 * (iMgradX * vy_local + iMgradY * vx_local);

            if (coordinates == SPHERICAL) {
                const auto iMM = _buffers.iMMDevice[dgi];
                e11.row(dgi) -= iMM * vy_local;
                e12.row(dgi) += 0.5 * iMM * vx_local;
            }
        });
}

// todo: should be KokkosMEVPDynamicsKernel::DeviceViewCG
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
void KokkosCGDynamicsKernel<DGadvection>::computeStressDivergence(const KokkosBuffers& _buffers,
    const KokkosMeshData& _meshData, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
{
    using CGVec = Eigen::Vector<Nextsim::FloatType, CGdof>;
    //    static PerfTimer timerDivZero("divZeroGPU");
    // static PerfTimer timerDivComp("divCompGPU");
    //   static PerfTimer timerDivDirichlet("divDirichletGPU");
    // zero buffers
    //   timerDivZero.start();
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, _buffers.dStressXDevice, 0.0);
    Kokkos::deep_copy(execSpace, _buffers.dStressYDevice, 0.0);
    //  execSpace.fence();
    //    timerDivZero.stop();

    //   timerDivComp.start();
    Kokkos::parallel_for(
        "computeStressDivergence", nx * ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            const DeviceIndex cx = idx % nx;
            const DeviceIndex cy = idx / nx;
            const DeviceIndex eid = cx + nx * cy;
            // only on ice!
            if (!_buffers.landMaskDevice.test(eid)) {
                return;
            }

            const auto s11 = makeEigenMap(_buffers.s11Device);
            const auto s12 = makeEigenMap(_buffers.s12Device);
            const auto s22 = makeEigenMap(_buffers.s22Device);

            const auto divS1 = _buffers.divS1Device[eid];
            const auto divS2 = _buffers.divS2Device[eid];
            CGVec tx = divS1 * s11.row(eid).transpose() + divS2 * s12.row(eid).transpose();
            CGVec ty = divS1 * s12.row(eid).transpose() + divS2 * s22.row(eid).transpose();

            if (coordinates == SPHERICAL) {
                const auto divM = _buffers.divMDevice[eid];
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
                    Kokkos::atomic_sub(&_buffers.dStressXDevice(dst_idx), tx(src_idx));
                    Kokkos::atomic_sub(&_buffers.dStressYDevice(dst_idx), ty(src_idx));
                }
            }
        });
    //    timerDivComp.stop();
    // set zero on the Dirichlet boundaries
    //   timerDivDirichlet.start();
    dirichletZero(_buffers.dStressXDevice, nx, ny, _meshData.dirichletDevice);
    dirichletZero(_buffers.dStressYDevice, nx, ny, _meshData.dirichletDevice);
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
    const DeviceViewCG& vDevice, const std::array<KokkosDeviceMapView, 4>& dirichlet,
    DeviceIndex nx, DeviceIndex ny)
{
    dirichletZero(uDevice, nx, ny, dirichletDevice);
    dirichletZero(vDevice, nx, ny, dirichletDevice);

    // TODO Periodic boundary conditions.
}

}