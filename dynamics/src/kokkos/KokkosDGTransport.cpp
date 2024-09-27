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

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const LocalVec<2> tangentTop = mesh.edgeVector(nx + mesh.nx + 1, nx + mesh.nx + 2);
            const LocalVec<EDGE_DOFS> vel2 = 0.5
                * (-tangentTop(0, 1) * topEdgeOfCell(velXDevice, cx)
                    + tangentTop(0, 0) * topEdgeOfCell(velYDevice, cx));
            addRowAtomic(normalVelXDevice, vel2, cx + mesh.nx);
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
template <int DG> void KokkosDGTransport<DG>::setTimeSteppingScheme(TimeSteppingScheme tss)
{
    timeSteppingScheme = tss;
}

/*************************************************************/
template <int DG>
void KokkosDGTransport<DG>::prepareAdvection(const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice)
{
    // todo: try interpolation in batches to fuse the kernels
    cG2DGInterpolator(velX, cgUDevice);
    cG2DGInterpolator(velY, cgVDevice);
    reinitNormalVelocityDevice(normalVelX, normalVelY, velX, velY, meshDevice);
}

/*************************************************************/
template <int DG> void KokkosDGTransport<DG>::step(FloatType dt, const DeviceViewDG& phiDevice)
{
/*    auto [phiHost, phiDevice] = makeKokkosDualView("phi", phi, true);
    Kokkos::deep_copy(phiDevice, phiHost);
    this->step_rk3(dt, phi);*/

    assert(phiDevice.size() == mesh.nx * mesh.ny * DG);

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

 /*   auto resGPU = phi;
    auto tmpRes1Host = makeKokkosHostView(resGPU);
    Kokkos::deep_copy(tmpRes1Host, phiDevice);
    compare("phi", phi, resGPU);*/
}

template <int DG>
void KokkosDGTransport<DG>::stepRK1(FloatType dt, const ConstDeviceViewDG& velXDevice,
    const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
    const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
    const DeviceViewDG& tmpRes1,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
    const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
    const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
    const KokkosMeshData& meshDevice)
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
    const KokkosMeshData& meshDevice)
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
                phiDevice(eid) += 0.5 * (tmpRes2(eid) - tmpRes1(eid));
            });
    } else {
        Kokkos::parallel_for(
            "add2", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                phiDevice(eid, ci) += 0.5 * (tmpRes2(eid, ci) - tmpRes1(eid, ci));
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
    const KokkosMeshData& meshDevice)
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
                tmpRes2(eid) = 0.25 * (tmpRes2(eid) + tmpRes1(eid)) + 0.75 * phiDevice(eid);
            });
    } else {
        Kokkos::parallel_for(
            "add2", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                tmpRes2(eid, ci)
                    = 0.25 * (tmpRes2(eid, ci) + tmpRes1(eid, ci)) + 0.75 * phiDevice(eid, ci);
            });
    }

    dGTransportOperatorDevice(dt, velX, velY, normalVelX, normalVelY, tmpRes2, tmpRes3,
        advectionCellTermXDevice, advectionCellTermYDevice, inverseDGMassMatrixDevice, meshDevice);

    if constexpr (DG == 1) {
        Kokkos::parallel_for(
            "add3", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex eid) {
                phiDevice(eid)
                    = 1.0 / 3.0 * phiDevice(eid) + 2.0 / 3.0 * (tmpRes2(eid) + tmpRes3(eid));
            });
    } else {
        Kokkos::parallel_for(
            "add3", phiDevice.size(), KOKKOS_LAMBDA(const DeviceIndex i) {
                const DeviceIndex eid = i / DG;
                const DeviceIndex ci = i % DG;
                phiDevice(eid, ci) = 1.0 / 3.0 * phiDevice(eid, ci)
                    + 2.0 / 3.0 * (tmpRes2(eid, ci) + tmpRes3(eid, ci));
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
                if (iy + 1 >= ny) {
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
                    dt * (std::max(vel, zero) * left + std::min(vel, zero) * right));
                Kokkos::atomic_add(&phiupDevice(c2),
                    dt * (std::max(vel, zero) * left + std::min(vel, zero) * right));
            });
    } else {
        const auto PSIe1D = PSIe<EDGE_DOFS, GP1D>;
        const auto PSIew1 = PSIe_w<DG, GP1D, 1>;
        const auto PSIew3 = PSIe_w<DG, GP1D, 3>;
        Kokkos::parallel_for(
            "edgeTermX", phiDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
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

                const LocalVec<GP1D> tmp = (velGauss.array().max(0)
                        * (rightEdgeOfCell(phiDevice, c1) * PSIe1D).array()
                    + velGauss.array().min(0) * (leftEdgeOfCell(phiDevice, c2) * PSIe1D).array());

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
    const KokkosMeshData& meshDevice)
{
    const auto PSIeE = PSIe<EDGE_DOFS, EDGE_DOFS>;
    // bot
    const auto PSIew0 = PSIe_w<DG, EDGE_DOFS, 0>;
    Kokkos::parallel_for(
        "dirichletTermBot", meshDevice.dirichletDevice[0].extent(0),
        KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = meshDevice.dirichletDevice[0][i];
            //    const DeviceIndex ix = eid % meshDevice.nx; // compute coordinates of element
            //    const DeviceIndex iy = eid / meshDevice.nx;
            const DeviceIndex c = eid;
            const DeviceIndex e = eid;

            const auto normalVelX = makeEigenMap(normalVelXDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelX.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp
                = (bottomEdgeOfCell(phiDevice, c) * PSIeE).array() * (-velGauss.array()).max(0);

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew0;
        });

    // right
    const auto PSIew1 = PSIe_w<DG, EDGE_DOFS, 1>;
    Kokkos::parallel_for(
        "dirichletTermRight", meshDevice.dirichletDevice[1].extent(0),
        KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = meshDevice.dirichletDevice[1][i];
            const DeviceIndex ix = eid % meshDevice.nx; // compute coordinates of element
            const DeviceIndex iy = eid / meshDevice.nx;
            const DeviceIndex c = eid;
            const DeviceIndex e = (meshDevice.nx + 1) * iy + ix + 1;

            const auto normalVelY = makeEigenMap(normalVelYDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelY.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp
                = (rightEdgeOfCell(phiDevice, c) * PSIeE).array() * velGauss.array().max(0);

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew1;
        });

    // top
    const auto PSIew2 = PSIe_w<DG, EDGE_DOFS, 2>;
    Kokkos::parallel_for(
        "dirichletTermTop", meshDevice.dirichletDevice[2].extent(0),
        KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = meshDevice.dirichletDevice[2][i];
            const DeviceIndex ix = eid % meshDevice.nx; // compute coordinates of element
            const DeviceIndex iy = eid / meshDevice.nx;
            const DeviceIndex c = eid;
            const DeviceIndex e = meshDevice.nx * (iy + 1) + ix;

            const auto normalVelX = makeEigenMap(normalVelXDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelX.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp
                = (topEdgeOfCell(phiDevice, c) * PSIeE).array() * velGauss.array().max(0);

            auto phiup = makeEigenMap(phiupDevice);
            phiup.row(c) -= dt * tmp * PSIew2;
        });

    // left
    const auto PSIew3 = PSIe_w<DG, EDGE_DOFS, 3>;
    Kokkos::parallel_for(
        "dirichletTermLeft", meshDevice.dirichletDevice[3].extent(0),
        KOKKOS_LAMBDA(const DeviceIndex i) {
            const DeviceIndex eid = meshDevice.dirichletDevice[3][i];
            const DeviceIndex ix = eid % meshDevice.nx; // compute coordinates of element
            const DeviceIndex iy = eid / meshDevice.nx;
            const DeviceIndex c = eid;
            const DeviceIndex e = (meshDevice.nx + 1) * iy + ix;

            const auto normalVelY = makeEigenMap(normalVelYDevice);
            const LocalVec<EDGE_DOFS> velGauss = normalVelY.row(e) * PSIeE;
            const LocalVec<EDGE_DOFS> tmp
                = (leftEdgeOfCell(phiDevice, c) * PSIeE).array() * (-velGauss.array()).max(0);

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
    const KokkosMeshData& meshDevice)
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