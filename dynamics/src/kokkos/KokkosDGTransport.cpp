/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosDGTransport.hpp"
#include "../include/ParametricTools.hpp"

namespace Nextsim {

template <int DG>
KokkosDGTransport<DG>::KokkosDGTransport(const ParametricMesh& smesh, const KokkosMesh& _meshDevice,
    const Interpolations::KokkosCG2DGInterpolator<DG, CGdegree>& _cG2DGInterpolator)
    : DGTransport<DG>(smesh)
    , meshDevice(_meshDevice)
    , cG2DGInterpolator(_cG2DGInterpolator)
    , timeSteppingScheme(TimeSteppingScheme::RK2)
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
    advectionCellTermXDevice = makeKokkosDeviceViewMap<MakeViewOptions::DEVICE_COPY>(
        "advectionCellTermX", this->parammap.AdvectionCellTermX);
    advectionCellTermYDevice = makeKokkosDeviceViewMap<MakeViewOptions::DEVICE_COPY>(
        "advectionCellTermY", this->parammap.AdvectionCellTermY);
    inverseDGMassMatrixDevice = makeKokkosDeviceViewMap<MakeViewOptions::DEVICE_COPY>(
        "inverseDGMassMatrix", this->parammap.InverseDGMassMatrix);
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

template <int Comps> using LocalVec = Eigen::Matrix<FloatType, 1, Comps>;

// dG0 (1 in cell, 1 on edge)

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
    return LocalVec<2>(cv(eid, 0) - FloatType(0.5) * cv(eid, 1), cv(eid, 2));
}

KOKKOS_IMPL_FUNCTION LocalVec<2> rightEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) + FloatType(0.5) * cv(eid, 1), cv(eid, 2));
}

KOKKOS_IMPL_FUNCTION LocalVec<2> bottomEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) - FloatType(0.5) * cv(eid, 2), cv(eid, 1));
}
KOKKOS_IMPL_FUNCTION LocalVec<2> topEdgeOfCell(const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return LocalVec<2>(cv(eid, 0) + FloatType(0.5) * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)

KOKKOS_IMPL_FUNCTION LocalVec<3> leftEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - FloatType(0.5) * cv(eid, 1) + FloatType(1. / 6.) * cv(eid, 3),
        cv(eid, 2) - FloatType(0.5) * cv(eid, 5), cv(eid, 4));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> rightEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + FloatType(0.5) * cv(eid, 1) + FloatType(1. / 6.) * cv(eid, 3),
        cv(eid, 2) + FloatType(0.5) * cv(eid, 5), cv(eid, 4));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> bottomEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - FloatType(0.5) * cv(eid, 2) + FloatType(1. / 6.) * cv(eid, 4),
        cv(eid, 1) - FloatType(0.5) * cv(eid, 5), cv(eid, 3));
}
KOKKOS_IMPL_FUNCTION LocalVec<3> topEdgeOfCell(const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + FloatType(0.5) * cv(eid, 2) + FloatType(1. / 6.) * cv(eid, 4),
        cv(eid, 1) + FloatType(0.5) * cv(eid, 5), cv(eid, 3));
}

// dG2+ (8 in cell, 3 on edge)

KOKKOS_IMPL_FUNCTION LocalVec<3> leftEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - FloatType(0.5) * cv(eid, 1) + FloatType(1. / 6.) * cv(eid, 3),
        cv(eid, 2) - FloatType(0.5) * cv(eid, 5) + FloatType(1. / 6.) * cv(eid, 6),
        cv(eid, 4) - FloatType(0.5) * cv(eid, 7));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> rightEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + FloatType(0.5) * cv(eid, 1) + FloatType(1. / 6.) * cv(eid, 3),
        cv(eid, 2) + FloatType(0.5) * cv(eid, 5) + FloatType(1. / 6.) * cv(eid, 6),
        cv(eid, 4) + FloatType(0.5) * cv(eid, 7));
}

KOKKOS_IMPL_FUNCTION LocalVec<3> bottomEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) - FloatType(0.5) * cv(eid, 2) + FloatType(1. / 6.) * cv(eid, 4),
        cv(eid, 1) - FloatType(0.5) * cv(eid, 5) + FloatType(1. / 6.) * cv(eid, 7),
        cv(eid, 3) - FloatType(0.5) * cv(eid, 6));
}
KOKKOS_IMPL_FUNCTION LocalVec<3> topEdgeOfCell(const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return LocalVec<3>(cv(eid, 0) + FloatType(0.5) * cv(eid, 2) + FloatType(1. / 6.) * cv(eid, 4),
        cv(eid, 1) + FloatType(0.5) * cv(eid, 5) + FloatType(1. / 6.) * cv(eid, 7),
        cv(eid, 3) + FloatType(0.5) * cv(eid, 6));
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
    const ConstDeviceViewDG& velYDevice, const KokkosMesh& mesh)
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

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const LocalVec<2> tangentLeft = mesh.edgeVector(ey, ey + mesh.nx + 1);
            const LocalVec<EDGE_DOFS> vel1 = FloatType(0.5)
                * (tangentLeft(0, 1) * leftEdgeOfCell(velXDevice, cy)
                    - tangentLeft(0, 0) * leftEdgeOfCell(velYDevice, cy));
            addRowAtomic(normalVelYDevice, vel1, ey);

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const LocalVec<2> tangentRight = mesh.edgeVector(ey + 1, ey + mesh.nx + 2);
            const LocalVec<EDGE_DOFS> vel2 = FloatType(0.5)
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
            const LocalVec<EDGE_DOFS> vel1 = FloatType(0.5)
                * (-tangentBottom(0, 1) * bottomEdgeOfCell(velXDevice, cx)
                    + tangentBottom(0, 0) * bottomEdgeOfCell(velYDevice, cx));
            addRowAtomic(normalVelXDevice, vel1, cx);

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const LocalVec<2> tangentTop = mesh.edgeVector(nx + mesh.nx + 1, nx + mesh.nx + 2);
            const LocalVec<EDGE_DOFS> vel2 = FloatType(0.5)
                * (-tangentTop(0, 1) * topEdgeOfCell(velXDevice, cx)
                    + tangentTop(0, 0) * topEdgeOfCell(velYDevice, cx));
            addRowAtomic(normalVelXDevice, vel2, cx + mesh.nx);
        });

    // TR 07.04.2025
    // correction of the normal velocity is needed on the mesh-boundary, where each edge
    // only has one neighbor. Above, we sum both sides, weighted with 0.5. Therefore,
    // we must scale the outer edges with 2. This is only used for inflow/outflow
    // Neumann boundaries.
    // bot
    Kokkos::parallel_for(
        "dirichletBot", mesh.nx, KOKKOS_LAMBDA(const DeviceIndex i) {
            auto normalVelX = makeEigenMap(normalVelXDevice);
            normalVelX.row(i) *= 2.0;
        });
    // right
    Kokkos::parallel_for(
        "dirichletRight", mesh.ny, KOKKOS_LAMBDA(const DeviceIndex i) {
            auto normalVelY = makeEigenMap(normalVelYDevice);
            normalVelY.row(i * (mesh.nx + 1) + mesh.nx) *= 2.0;
        });
    // top
    Kokkos::parallel_for(
        "dirichletTop", mesh.nx, KOKKOS_LAMBDA(const DeviceIndex i) {
            auto normalVelX = makeEigenMap(normalVelXDevice);
            normalVelX.row(i + mesh.nx * mesh.ny) *= 2.0;
        });
    // Left
    Kokkos::parallel_for(
        "dirichletLeft", mesh.ny, KOKKOS_LAMBDA(const DeviceIndex i) {
            auto normalVelY = makeEigenMap(normalVelYDevice);
            normalVelY.row(i * (mesh.nx + 1)) *= 2.0;
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
template <int DG> void KokkosDGTransport<DG>::setTimeSteppingScheme(TimeSteppingScheme tss)
{
    timeSteppingScheme = tss;
}

/*************************************************************/
template <int DG>
void KokkosDGTransport<DG>::prepareAdvection(
    const ConstKokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const ConstKokkosDeviceView<CGVector<CGdegree>>& cgVDevice)
{
    // todo: try interpolation in batches to fuse the kernels
    cG2DGInterpolator(velX, cgUDevice);
    cG2DGInterpolator(velY, cgVDevice);
    reinitNormalVelocityDevice(normalVelX, normalVelY, velX, velY, meshDevice);
}

/*************************************************************/
template <int DG> void KokkosDGTransport<DG>::step(FloatType dt, const DeviceViewDG& phiDevice)
{
    assert(phiDevice.size() == meshDevice.nx * meshDevice.ny * DG);

    switch (timeSteppingScheme) {
    case TimeSteppingScheme::RK1:
        stepRK1(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1,
            advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice,
            meshDevice);
        break;
    case TimeSteppingScheme::RK2:
        stepRK2(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1, tmpRes2,
            advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice,
            meshDevice);
        break;
    case TimeSteppingScheme::RK3:
        stepRK3(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1, tmpRes2, tmpRes3,
            advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice,
            meshDevice);
        break;
    }
}

template <int DG>
void KokkosDGTransport<DG>::stepRK1(FloatType dt, const ConstDeviceViewDG& velXDevice,
    const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
    const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
    const DeviceViewDG& tmpRes1,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
    const KokkosMesh& meshDevice)
{
    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add", phiDevice.size(),
            KOKKOS_LAMBDA(const DeviceIndex eid) { phiDevice(eid) += tmpRes1(eid); });
    } else {
        Kokkos::parallel_for(
            "add", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                phiDevice(eid, ci) += tmpRes1(eid, ci);
            });
    }
}

template <int DG>
void KokkosDGTransport<DG>::stepRK2(FloatType dt, const ConstDeviceViewDG& velXDevice,
    const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
    const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
    const DeviceViewDG& tmpRes1, const DeviceViewDG& tmpRes2,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
    const KokkosMesh& meshDevice)
{
    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add", phiDevice.size(),
            KOKKOS_LAMBDA(const DeviceIndex eid) { phiDevice(eid) += tmpRes1(eid); });
    } else {
        Kokkos::parallel_for(
            "add", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                phiDevice(eid, ci) += tmpRes1(eid, ci);
            });
    }

    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes2,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add2", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex eid) {
                phiDevice(eid) += FloatType(0.5) * (tmpRes2(eid) - tmpRes1(eid));
            });
    } else {
        Kokkos::parallel_for(
            "add2", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                phiDevice(eid, ci) += FloatType(0.5) * (tmpRes2(eid, ci) - tmpRes1(eid, ci));
            });
    }
}

template <int DG>
void KokkosDGTransport<DG>::stepRK3(FloatType dt, const ConstDeviceViewDG& velXDevice,
    const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
    const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
    const DeviceViewDG& tmpRes1, const DeviceViewDG& tmpRes2, const DeviceViewDG& tmpRes3,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
    const KokkosMesh& meshDevice)
{
    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, phiDevice, tmpRes1,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add", phiDevice.size(),
            KOKKOS_LAMBDA(const DeviceIndex eid) { tmpRes1(eid) += phiDevice(eid); });
    } else {
        Kokkos::parallel_for(
            "add", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                tmpRes1(eid, ci) += phiDevice(eid, ci);
            });
    }

    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, tmpRes1, tmpRes2,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add2", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex eid) {
                tmpRes2(eid)
                    = FloatType(0.25) * (tmpRes2(eid) + tmpRes1(eid)) + 0.75 * phiDevice(eid);
            });
    } else {
        Kokkos::parallel_for(
            "add2", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                tmpRes2(eid, ci) = FloatType(0.25) * (tmpRes2(eid, ci) + tmpRes1(eid, ci))
                    + FloatType(0.75) * phiDevice(eid, ci);
            });
    }

    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, tmpRes2, tmpRes3,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add3", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex eid) {
                phiDevice(eid) = FloatType(1.0 / 3.0) * phiDevice(eid)
                    + FloatType(2.0 / 3.0) * (tmpRes2(eid) + tmpRes3(eid));
            });
    } else {
        Kokkos::parallel_for(
            "add3", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                phiDevice(eid, ci) = FloatType(1.0 / 3.0) * phiDevice(eid, ci)
                    + FloatType(2.0 / 3.0) * (tmpRes2(eid, ci) + tmpRes3(eid, ci));
            });
    }
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

/*************************************************************/
template <int DG>
void KokkosDGTransport<DG>::addEdgeXTermsDevice(FloatType dt,
    const ConstDeviceViewEdge& normalVelXDevice, const ConstDeviceViewDG& phiDevice,
    const DeviceViewDG& phiupDevice, const ConstDeviceBitset& landMaskDevice, DeviceIndex nx,
    DeviceIndex ny)
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
                    dt * (Kokkos::fmax(vel, zero) * bottom + Kokkos::fmin(vel, zero) * top));
                Kokkos::atomic_add(&phiupDevice(c2),
                    dt * (Kokkos::fmax(vel, zero) * bottom + Kokkos::fmin(vel, zero) * top));
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
                if (iy + 1 >= ny) {
                    return;
                }

                const auto normalVelX = makeEigenMap(normalVelXDevice);
                const LocalVec<GP1D> velGauss = normalVelX.row(ie) * PSIe1D;

                const LocalVec<GP1D> tmp = (velGauss.array().max(FloatType(0))
                        * (topEdgeOfCell(phiDevice, c1) * PSIe1D).array()
                    + velGauss.array().min(FloatType(0))
                        * (bottomEdgeOfCell(phiDevice, c2) * PSIe1D).array());

                addRowAtomic(phiupDevice, (-dt * tmp * PSIew2).eval(), c1);
                addRowAtomic(phiupDevice, (dt * tmp * PSIew0).eval(), c2);
            });
    }
}

template <int DG>
void KokkosDGTransport<DG>::addEdgeYTermsDevice(FloatType dt,
    const ConstDeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& phiDevice,
    const DeviceViewDG& phiupDevice, const ConstDeviceBitset& landMaskDevice, DeviceIndex nx,
    DeviceIndex ny)
{
    // branch needs to be outside because it interferes with lambda implicit captures
    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "edgeTermY", phiDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
                const DeviceIndex ix = eid % nx;
                const DeviceIndex iy = eid / nx;
                const DeviceIndex c1 = eid;
                const DeviceIndex c2 = c1 + 1;
                // first index of inner velocity in row
                const DeviceIndex ie = iy * (nx + 1) + 1 + ix;

                if (!landMaskDevice.test(c1) || !landMaskDevice.test(c2)) {
                    return;
                }
                // only inner edges
                if (iy >= ny) {
                    return;
                }

                const FloatType left = phiDevice(c1);
                const FloatType right = phiDevice(c2);
                const FloatType vel = normalVelYDevice(ie);

                // max and min would not compile if the float literal type is not the same as
                // FloatType
                constexpr FloatType zero = 0;
                Kokkos::atomic_sub(&phiupDevice(c1),
                    dt * (Kokkos::fmax(vel, zero) * left + Kokkos::fmin(vel, zero) * right));
                Kokkos::atomic_add(&phiupDevice(c2),
                    dt * (Kokkos::fmax(vel, zero) * left + Kokkos::fmin(vel, zero) * right));
            });
    } else {
        const auto PSIe1D = PSIe<EDGE_DOFS, GP1D>;
        const auto PSIew1 = PSIe_w<DG, GP1D, 1>;
        const auto PSIew3 = PSIe_w<DG, GP1D, 3>;
        Kokkos::parallel_for(
            "edgeTermY", phiDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
                const DeviceIndex ix = eid % nx;
                const DeviceIndex iy = eid / nx;
                const DeviceIndex c1 = eid;
                const DeviceIndex c2 = c1 + 1;
                // first index of inner velocity in row
                const DeviceIndex ie = iy * (nx + 1) + 1 + ix;

                if (!landMaskDevice.test(c1) || !landMaskDevice.test(c2)) {
                    return;
                }
                // only inner edges
                if (ix + 1 >= nx) {
                    return;
                }

                const auto normalVelY = makeEigenMap(normalVelYDevice);
                const LocalVec<GP1D> velGauss = normalVelY.row(ie) * PSIe1D;

                const LocalVec<GP1D> tmp = (velGauss.array().max(FloatType(0))
                        * (rightEdgeOfCell(phiDevice, c1) * PSIe1D).array()
                    + velGauss.array().min(FloatType(0))
                        * (leftEdgeOfCell(phiDevice, c2) * PSIe1D).array());

                addRowAtomic(phiupDevice, (-dt * tmp * PSIew1).eval(), c1);
                addRowAtomic(phiupDevice, (dt * tmp * PSIew3).eval(), c2);
            });
    }
}

/*************************************************************/
template <int DG>
void KokkosDGTransport<DG>::addBoundaryTermsDevice(FloatType dt,
    const ConstDeviceViewEdge& normalVelXDevice, const ConstDeviceViewEdge& normalVelYDevice,
    const ConstDeviceViewDG& phiDevice, const DeviceViewDG& phiupDevice,
    const KokkosMesh& meshDevice)
{
    const auto PSIeE = PSIe<EDGE_DOFS, EDGE_DOFS>;
    // bot
    const auto PSIew0 = PSIe_w<DG, EDGE_DOFS, 0>;
    Kokkos::parallel_for(
        "dirichletTermBot", meshDevice.nx, KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex c = i;
            if (!meshDevice.landMaskDevice.test(c)) {
                return;
            }
            const DeviceIndex e = i;

            const auto normalVelX = makeEigenMap(normalVelXDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelX.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp = (bottomEdgeOfCell(phiDevice, c) * PSIeE).array()
                * (-velGauss.array()).max(FloatType(0));

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew0;
        });

    // right
    const auto PSIew1 = PSIe_w<DG, EDGE_DOFS, 1>;
    Kokkos::parallel_for(
        "dirichletTermRight", meshDevice.ny, KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex c = (i + 1) * meshDevice.nx - 1;
            if (!meshDevice.landMaskDevice.test(c)) {
                return;
            }
            const DeviceIndex e = i * (meshDevice.nx + 1) + meshDevice.nx;

            const auto normalVelY = makeEigenMap(normalVelYDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelY.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp = (rightEdgeOfCell(phiDevice, c) * PSIeE).array()
                * velGauss.array().max(FloatType(0));

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew1;
        });

    // top
    const auto PSIew2 = PSIe_w<DG, EDGE_DOFS, 2>;
    Kokkos::parallel_for(
        "dirichletTermTop", meshDevice.nx, KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex c = i + meshDevice.nx * (meshDevice.ny - 1);
            if (!meshDevice.landMaskDevice.test(c)) {
                return;
            }
            const DeviceIndex e = i + meshDevice.nx * meshDevice.ny;

            const auto normalVelX = makeEigenMap(normalVelXDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelX.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp = (topEdgeOfCell(phiDevice, c) * PSIeE).array()
                * velGauss.array().max(FloatType(0));

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew2;
        });

    // left
    const auto PSIew3 = PSIe_w<DG, EDGE_DOFS, 3>;
    Kokkos::parallel_for(
        "dirichletTermLeft", meshDevice.ny, KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex c = i * meshDevice.nx;
            if (!meshDevice.landMaskDevice.test(c)) {
                return;
            }
            const DeviceIndex e = i * (meshDevice.nx + 1);

            const auto normalVelY = makeEigenMap(normalVelYDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelY.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp = (leftEdgeOfCell(phiDevice, c) * PSIeE).array()
                * (-velGauss.array()).max(FloatType(0));

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew3;
        });
}

template <int DG>
void KokkosDGTransport<DG>::dGTransportOperatorDevice(FloatType dt,
    const ConstDeviceViewDG& velXDevice, const ConstDeviceViewDG& velYDevice,
    const ConstDeviceViewEdge& normalVelXDevice, const ConstDeviceViewEdge& normalVelYDevice,
    const ConstDeviceViewDG& phiDevice, const DeviceViewDG& phiupDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
    const KokkosMesh& meshDevice)
{
    Kokkos::deep_copy(phiupDevice, 0.0);

    addCellTermsDevice(dt, velXDevice, velYDevice, phiDevice, phiupDevice, advectionCellTermXDevice,
        advectionCellTermYDevice, meshDevice.landMaskDevice);

    addEdgeXTermsDevice(dt, normalVelXDevice, phiDevice, phiupDevice, meshDevice.landMaskDevice,
        meshDevice.nx, meshDevice.ny);

    addEdgeYTermsDevice(dt, normalVelYDevice, phiDevice, phiupDevice, meshDevice.landMaskDevice,
        meshDevice.nx, meshDevice.ny);

    addBoundaryTermsDevice(
        dt, normalVelXDevice, normalVelYDevice, phiDevice, phiupDevice, meshDevice);

    Kokkos::parallel_for(
        "inverseMap", meshDevice.nx * meshDevice.ny, KOKKOS_LAMBDA(const DeviceIndex eid) {
            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(eid) = inverseDGMassMatrixDevice[eid] * phiup.row(eid).transpose();
        });
}

template class KokkosDGTransport<1>;
template class KokkosDGTransport<3>;
template class KokkosDGTransport<6>;
template class KokkosDGTransport<8>;
}
