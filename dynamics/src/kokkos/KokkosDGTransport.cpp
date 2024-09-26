/*!
 * @file KokkosDGTransport.cpp
 * @date September 12, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosDGTransport.hpp"
#include "../include/ParametricTools.hpp"

namespace Nextsim {

template <int DG>
KokkosDGTransport<DG>::KokkosDGTransport(const ParametricMesh& smesh,
    const KokkosMeshData& _meshDevice,
    const Interpolations::KokkosCG2DGInterpolator<DG, CGdegree>& _cG2DGInterpolator)
    : DGTransport<DG>(smesh)
    , meshDevice(_meshDevice)
    , cG2DGInterpolator(_cG2DGInterpolator)
{
    // todo: initialize without hostTransport
    velX = makeKokkosDeviceView("velX", this->velx);
    velY = makeKokkosDeviceView("velY", this->vely);

    normalVelX = makeKokkosDeviceView("normalVelX", this->normalvel_X);
    normalVelY = makeKokkosDeviceView("normalVelY", this->normalvel_Y);

    tmpRes1 = makeKokkosDeviceView("tmpRes1", this->tmp1);
    tmpRes2 = makeKokkosDeviceView("tmpRes2", this->tmp2);
    tmpRes3 = makeKokkosDeviceView("tmpRes3", this->tmp3);

    // parametric map
    advectionCellTermXDevice
        = makeKokkosDeviceViewMap("advectionCellTermX", this->parammap.AdvectionCellTermX, true);
    advectionCellTermYDevice
        = makeKokkosDeviceViewMap("advectionCellTermY", this->parammap.AdvectionCellTermY, true);

    inverseDGMassMatrixDevice
        = makeKokkosDeviceViewMap("inverseDGMassMatrix", this->parammap.InverseDGMassMatrix, true);
}

//! returns the localization of the cell vector to the edges
/*!
 * writes the cell-basis on the edges in the edge basis
 *
 * CELL:
 * DGdegree 0:      1
 * DGdegree 1-2:    + (x-1/2), (y-1/2),
 * DGdegree 3-5:    + (x-1/2)^2-1/12, (y-1/2)^2-1/12, (x-1/2)(y-1/2)
 * DGdegree 6-7:    + (y-1/2)(x-1/2)^2-1/12, (x-1/2)(y-1/2)^2-1/12
 *
 * EDGE:
 * DGdegree 0:    1
 * DGdegree 1: +  (t-1/2)
 * DGdegree 2: +  (t-1/2)^2-1/12
 *
 */

template <int DG> using ConstDeviceViewDG = ConstKokkosDeviceView<DGVector<DG>>;
/*
template <int DG>
LocalVec<EDGE_DOFS> leftEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);
template <int DG>
LocalVec<EDGE_DOFS> rightEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);
template <int DG>
LocalVec<EDGE_DOFS> bottomEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);
template <int DG>
LocalVec<EDGE_DOFS> topEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);*/
/*
namespace Details {
    // extract the DGDegree from an Eigen or Kokkos array
    template <typename DGVec> struct DGDegree;

    // Kokkos
    template <typename T, int DG, class... Properties>
    struct DGDegree<Kokkos::View<T* [DG], Properties...>> {
        static constexpr int value = DG;
    };

    // Eigen
    template <typename T, int DG, int Options>
    struct DGDegree<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, DG, Options>> {
        static constexpr int value = DG;
    };
}

template <typename DGVec>
static KOKKOS_IMPL_FUNCTION auto leftEdgeOfCell(const DGVec& cv, DeviceIndex eid)
{
    constexpr int DG = Details::DGDegree<DGVec>::value;
    if constexpr (DG == 1) {
        return LocalVec<1>(cv(eid, 0));
    } else if constexpr (DG == 3){
        return LocalVec<2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
    } else if constexpr (DG == 6) {
    }
    } else {
        return LocalVec<1>(cv(eid, 0));
    }
}*/

// dG0 (1 in cell, 1 on edge)

template <int Comps> using LocalVec = Eigen::Matrix<FloatType, 1, Comps>;

KOKKOS_IMPL_FUNCTION LocalVec<1> leftEdgeOfCell(const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return LocalVec<1>(cv(eid));
}

KOKKOS_IMPL_FUNCTION LocalVec<1> rightEdgeOfCell(const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return LocalVec<1>(cv(eid));
}

KOKKOS_IMPL_FUNCTION LocalVec<1> bottomEdgeOfCell(const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return LocalVec<1>(cv(eid));
}
KOKKOS_IMPL_FUNCTION LocalVec<1> topEdgeOfCell(const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return LocalVec<1>(cv(eid));
}

// dG1 (3 in cell, 2 on edge)

KOKKOS_IMPL_FUNCTION LocalVec<2> leftEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}

KOKKOS_IMPL_FUNCTION LocalVec<2> rightEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}

KOKKOS_IMPL_FUNCTION LocalVec<2> bottomEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
KOKKOS_IMPL_FUNCTION LocalVec<2> topEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)

KOKKOS_IMPL_FUNCTION LocalVec<3> leftEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5), cv(eid, 4));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> rightEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5), cv(eid, 4));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> bottomEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5), cv(eid, 3));
}
KOKKOS_IMPL_FUNCTION LocalVec<3> topEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5), cv(eid, 3));
}

// dG2+ (8 in cell, 3 on edge)

KOKKOS_IMPL_FUNCTION LocalVec<3> leftEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) - 0.5 * cv(eid, 7));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> rightEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) + 0.5 * cv(eid, 7));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> bottomEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) - 0.5 * cv(eid, 6));
}
KOKKOS_IMPL_FUNCTION LocalVec<3> topEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) + 0.5 * cv(eid, 6));
}

// add local vector to global row atomically
template <typename T, int Comps, int Options, class... Properties>
KOKKOS_IMPL_FUNCTION static void addRowAtomic(const Kokkos::View<T* [Comps], Properties...>& dst,
    const Eigen::Matrix<T, 1, Comps, Options>& src, DeviceIndex row)
{
    for (int i = 0; i < Comps; ++i) {
        Kokkos::atomic_add(&dst(row, i), src(i));
    }
}
// special case because a 1D kokkos view expects only one index
template <typename T, int Options, class... Properties>
KOKKOS_IMPL_FUNCTION static void addRowAtomic(const Kokkos::View<T*, Properties...>& dst,
    const Eigen::Matrix<T, 1, 1, Options>& src, DeviceIndex row)
{
    Kokkos::atomic_add(&dst(row), src(0));
}

template <int DG>
void KokkosDGTransport<DG>::reinitNormalVelocityDevice(const DeviceViewEdge& normalVelXDevice,
    const DeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& velXDevice,
    const ConstDeviceViewDG& velYDevice, const KokkosMeshData& mesh)
{
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, normalVelXDevice, 0.0);
    Kokkos::deep_copy(execSpace, normalVelYDevice, 0.0);

    // average the velocity to the Y-edges
    Kokkos::parallel_for(
        "averageVelocityY", mesh.nx * mesh.ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            if (!mesh.landMaskDevice.test(idx)) {
                return;
            }

            const DeviceIndex ix = idx % mesh.nx;
            const DeviceIndex iy = idx / mesh.nx;
            //   |     |
            // --*-----*--
            //  ey  cy |
            //   |     |
            // -ey-----*--
            //   |     |
            // first edge-index and node-index in row
            const DeviceIndex ey = ix + iy * (mesh.nx + 1);
            const DeviceIndex cy = ix + iy * mesh.nx; // first cell index in row

            // todo: check if this write has a race-condition
            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const LocalVec<2> tangentLeft = mesh.edgeVector(ey, ey + mesh.nx + 1);
            const LocalVec<EDGE_DOFS> vel1 = 0.5
                * (tangentLeft(0, 1) * leftEdgeOfCell(velXDevice, cy)
                    - tangentLeft(0, 0) * leftEdgeOfCell(velYDevice, cy));
            addRowAtomic(normalVelYDevice, vel1, ey);

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const LocalVec<2> tangentRight = mesh.edgeVector(ey + 1, ey + mesh.nx + 2);
            const LocalVec<EDGE_DOFS> vel2 = 0.5
                * (tangentRight(0, 1) * rightEdgeOfCell(velXDevice, cy)
                    - tangentRight(0, 0) * rightEdgeOfCell(velYDevice, cy));
            addRowAtomic(normalVelYDevice, vel2, ey + 1);
            // we need an adjustment along the boundaries.. This is done later on.
        });

    // average the velocity to the X-edges
    Kokkos::parallel_for(
        "averageVelocityX", mesh.nx * mesh.ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            if (!mesh.landMaskDevice.test(idx)) {
                return;
            }
            const DeviceIndex ix = idx % mesh.nx;
            const DeviceIndex iy = idx / mesh.nx;
            //   |     |
            // --*-----*--
            //   |  cx |
            //   |     |
            // -mesh.nx-ex--*--
            //   |     |

            const DeviceIndex cx = ix + iy * mesh.nx; // first edge-index and cell-index
            const DeviceIndex nx = ix + iy * (mesh.nx + 1); // first cell index in row

            // un-normed tangent vector of bottom edge (pointing right). normal is (-y,x)
            const LocalVec<2> tangentBottom = mesh.edgeVector(nx, nx + 1);
            const LocalVec<EDGE_DOFS> vel1 = 0.5
                * (-tangentBottom(0, 1) * bottomEdgeOfCell(velXDevice, cx)
                    + tangentBottom(0, 0) * bottomEdgeOfCell(velYDevice, cx));
            addRowAtomic(normalVelXDevice, vel1, cx);

            /*    normalVelX.row(cx) += 0.5
                    * (-tangentBottom(0, 1) * bottomEdgeOfCell(velXDevice, cx)
                        + tangentBottom(0, 0) * bottomEdgeOfCell(velYDevice, cx));*/

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const LocalVec<2> tangentTop = mesh.edgeVector(nx + mesh.nx + 1, nx + mesh.nx + 2);
            const LocalVec<EDGE_DOFS> vel2 = 0.5
                * (-tangentTop(0, 1) * topEdgeOfCell(velXDevice, cx)
                    + tangentTop(0, 0) * topEdgeOfCell(velYDevice, cx));
            addRowAtomic(normalVelXDevice, vel2, cx + mesh.nx);

            /*   normalVelX.row(cx + mesh.nx) += 0.5
                   * (-tangentTop(0, 1) * topEdgeOfCell(velXDevice, cx)
                       + tangentTop(0, 0) * topEdgeOfCell(velYDevice, cx));*/
        });

    // Take care of the boundaries. Usually, the normal velocity is the average velocity
    // from left and from the right. Hence, we get the factor 0.5 above. At boundaries,
    // the normal is set only once, from the inside. These edges must be scaled with 2.0
    // todo: measure if capturing the whole KokkosMesh makes a difference here
    // bot
    Kokkos::parallel_for(
        "dirichletBot", mesh.dirichletDevice[0].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = mesh.dirichletDevice[0][i];
            const DeviceIndex ix = eid % mesh.nx; // compute coordinates of element
            const DeviceIndex iy = eid / mesh.nx;
            // todo: check if working directly on normalVelXDevice is faster
            auto normalVelX = makeEigenMap(normalVelXDevice);
            normalVelX.row(mesh.nx * iy + ix) *= 2.0;
        });
    // right
    Kokkos::parallel_for(
        "dirichletRight", mesh.dirichletDevice[1].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = mesh.dirichletDevice[1][i];
            const DeviceIndex ix = eid % mesh.nx; // compute coordinates of element
            const DeviceIndex iy = eid / mesh.nx;
            auto normalVelY = makeEigenMap(normalVelYDevice);
            normalVelY.row((mesh.nx + 1) * iy + ix + 1) *= 2.0;
        });
    // top
    Kokkos::parallel_for(
        "dirichletTop", mesh.dirichletDevice[2].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = mesh.dirichletDevice[2][i];
            const DeviceIndex ix = eid % mesh.nx; // compute coordinates of element
            const DeviceIndex iy = eid / mesh.nx;
            auto normalVelX = makeEigenMap(normalVelXDevice);
            normalVelX.row(mesh.nx * (iy + 1) + ix) *= 2.0;
        });
    // left
    Kokkos::parallel_for(
        "dirichletLeft", mesh.dirichletDevice[3].extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = mesh.dirichletDevice[3][i];
            const DeviceIndex ix = eid % mesh.nx; // compute coordinates of element
            const DeviceIndex iy = eid / mesh.nx;
            auto normalVelY = makeEigenMap(normalVelYDevice);
            normalVelY.row((mesh.nx + 1) * iy + ix) *= 2.0;
        });
}

template <typename Mat> void compare(const std::string& name, const Mat& m1, const Mat& m2)
{
    FloatType normRef = m1.norm();
    FloatType normDiff = (m1 - m2).norm();
    std::cout << name << " - abs: " << normDiff << ", rel: " << normDiff / normRef
              << ", norm: " << normRef << std::endl;
}

/*************************************************************/
template <int DG>
// template <int CG>
void KokkosDGTransport<DG>::prepareAdvection(const CGVector<CGdegree>& cgU,
    const CGVector<CGdegree>& cgV, const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice)
{
    DGTransport<DG>::prepareAdvection(cgU, cgV);
    /*    auto velXHost = makeKokkosHostView(this->velx);
        auto velYHost = makeKokkosHostView(this->vely);
        Kokkos::deep_copy(velX, velXHost);
        Kokkos::deep_copy(velY, velYHost);*/
    // todo: try interpolation in batches to fuse the kernels
    cG2DGInterpolator(velX, cgUDevice);
    cG2DGInterpolator(velY, cgVDevice);
    reinitNormalVelocityDevice(normalVelX, normalVelY, velX, velY, meshDevice);

    auto tempX = this->normalvel_X;
    auto tempY = this->normalvel_Y;
    auto normalVelXHost = makeKokkosHostView(this->normalvel_X);
    auto normalVelYHost = makeKokkosHostView(this->normalvel_Y);
    Kokkos::deep_copy(normalVelXHost, normalVelX);
    Kokkos::deep_copy(normalVelYHost, normalVelY);
    compare("normalVelX", tempX, this->normalvel_X);
    compare("normalVelY", tempY, this->normalvel_Y);
}

/*************************************************************/
template <int DG> void KokkosDGTransport<DG>::step(FloatType dt, DGVector<DG>& phi)
{
    this->step_rk1(dt, phi);

    auto [phiHost, phiDevice] = makeKokkosDualView("phi", phi, true);
    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);
    auto resGPU = this->tmp1;
    auto tmpRes1Host = makeKokkosHostView(resGPU);
    Kokkos::deep_copy(tmpRes1Host, tmpRes1);
    compare("tmp1", this->tmp1, resGPU);
}

/*************************************************************/
template <int DG>
void KokkosDGTransport<DG>::addCellTermsDevice(FloatType dt, const ConstDeviceViewDG& velXDevice,
    const ConstDeviceViewDG& velYDevice, const ConstDeviceViewDG& phiDevice,
    const DeviceViewDG& phiupDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const ConstDeviceBitset& landMaskDevice)
{
    // allow for explicit capture of the constant data
    // In principal constant memory is used to store the functor / args of the kernel.
    // However, it means that the data needs to be copied for every kernel invocation and there
    // is a likely a performance penalty for large functors.
    // Kokkos uses shared memory instead of constant memory to store the functor if >512 bytes
    // (https://github.com/kokkos/kokkos/issues/606) which could be the case if just a single PSI
    // matrix is captured. todo: measure which way is faster
    const auto PSIDG = PSI<DG, GP1D>;
    if constexpr (DG > 1) {
        Kokkos::parallel_for(
            "cellTerm", phiDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
                if (!landMaskDevice.test(eid)) {
                    return;
                }

                const auto vx = makeEigenMap(velXDevice);
                const auto vy = makeEigenMap(velYDevice);
                const auto phi = makeEigenMap(phiDevice);

                //!< velocity in GP
                const LocalVec<GP> vxGauss = vx.row(eid) * PSIDG;
                const LocalVec<GP> vyGauss = vy.row(eid) * PSIDG;
                const LocalVec<GP> phiGauss = (phi.row(eid) * PSIDG).array();

                const auto& advectionCellTermX = advectionCellTermXDevice[eid];
                const auto& advectionCellTermY = advectionCellTermYDevice[eid];
                auto phiup = makeEigenMap(phiupDevice);
                phiup.row(eid) += dt
                    * (advectionCellTermX.array().rowwise() * vxGauss.array()
                        + advectionCellTermY.array().rowwise() * vyGauss.array())
                          .matrix()
                    * phiGauss.transpose();
            });
    }
}

template <int DG>
void KokkosDGTransport<DG>::addEdgeXTermsDevice(FloatType dt,
    const ConstDeviceViewEdge& normalVelXDevice, const ConstDeviceViewEdge& normalVelYDevice,
    const ConstDeviceViewDG& phiDevice, const DeviceViewDG& phiupDevice,
    const ConstDeviceBitset& landMaskDevice, DeviceIndex nx, DeviceIndex ny)
{
    // branch needs to be outside because it interferes with lambda implicit captures
    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "edgeTermX", phiDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
                const DeviceIndex ix = eid % nx;
                const DeviceIndex iy = eid / nx;
                const DeviceIndex c1 = eid;
                const DeviceIndex c2 = c1 + nx;
                const DeviceIndex ie = ix + nx + ny * nx;

                if (!landMaskDevice.test(c1) || !landMaskDevice.test(c2)) {
                    return;
                }
                // only inner edges
                if (iy >= ny) {
                    return;
                }

                const FloatType bottom = phiDevice(c1);
                const FloatType top = phiDevice(c2);
                const FloatType vel = normalVelXDevice(ie);

                // max and min would not compile if the float literal type is not the same as
                // FloatType
                constexpr FloatType zero = 0;
                Kokkos::atomic_sub(&phiupDevice(c1),
                    dt * (std::max(vel, zero) * bottom + std::min(vel, zero) * top));
                Kokkos::atomic_add(&phiupDevice(c2),
                    dt * (std::max(vel, zero) * bottom + std::min(vel, zero) * top));
            });
    } else {
        const auto PSIe1D = PSIe<EDGE_DOFS, GP1D>;
        const auto PSIew2 = PSIe_w<DG, GP1D, 2>;
        const auto PSIew0 = PSIe_w<DG, GP1D, 0>;
        Kokkos::parallel_for(
            "edgeTermX", phiDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
                const DeviceIndex ix = eid % nx;
                const DeviceIndex iy = eid / nx;
                const DeviceIndex c1 = eid;
                const DeviceIndex c2 = c1 + nx;
                const DeviceIndex ie = ix + nx + iy * nx;

                if (!landMaskDevice.test(c1) || !landMaskDevice.test(c2)) {
                    return;
                }
                // only inner edges
                if (iy+1 >= ny) {
                    return;
                }

                const auto normalVelX = makeEigenMap(normalVelXDevice);
                const LocalVec<GP1D> velGauss = normalVelX.row(ie) * PSIe1D;

                const LocalVec<GP1D> tmp = (velGauss.array().max(0)
                        * (topEdgeOfCell(phiDevice, c1) * PSIe1D).array()
                    + velGauss.array().min(0) * (bottomEdgeOfCell(phiDevice, c2) * PSIe1D).array());

                addRowAtomic(phiupDevice, (-dt * tmp * PSIew2).eval(), c1);
                addRowAtomic(phiupDevice, (dt * tmp * PSIew0).eval(), c2);
            });
    }
}

template <int DG>
void KokkosDGTransport<DG>::dGTransportOperatorDevice(FloatType dt,
    const ConstDeviceViewDG& velXDevice, const ConstDeviceViewDG& velYDevice,
    const ConstDeviceViewEdge& normalVelXDevice, const ConstDeviceViewEdge& normalVelYDevice,
    const ConstDeviceViewDG& phiDevice, const DeviceViewDG& phiupDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
    const KokkosMeshData& meshDevice)
{
    Kokkos::deep_copy(phiupDevice, 0.0);

    addCellTermsDevice(dt, velXDevice, velYDevice, phiDevice, phiupDevice, advectionCellTermXDevice,
        advectionCellTermYDevice, meshDevice.landMaskDevice);

    addEdgeXTermsDevice(dt, normalVelXDevice, normalVelYDevice, phiDevice, phiupDevice,
        meshDevice.landMaskDevice, meshDevice.nx, meshDevice.ny);
}

template class KokkosDGTransport<1>;
template class KokkosDGTransport<3>;
template class KokkosDGTransport<6>;
template class KokkosDGTransport<8>;
/*
template void KokkosDGTransport<1>::prepareAdvection(const CGVector<CGdegree>& cgU,
    const CGVector<CGdegree>& cgV, const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice);
template void KokkosDGTransport<3>::prepareAdvection(const CGVector<CGdegree>& cgU,
    const CGVector<CGdegree>& cgV, const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice);
template void KokkosDGTransport<6>::prepareAdvection(const CGVector<CGdegree>& cgU,
    const CGVector<CGdegree>& cgV, const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice);
template void KokkosDGTransport<8>::prepareAdvection(const CGVector<CGdegree>& cgU,
    const CGVector<CGdegree>& cgV, const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice);*/
}