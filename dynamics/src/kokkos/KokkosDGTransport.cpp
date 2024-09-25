/*!
 * @file KokkosDGTransport.cpp
 * @date September 12, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosDGTransport.hpp"
#include "../include/ParametricTools.hpp"

namespace Nextsim {

template <int DG>
KokkosDGTransport<DG>::KokkosDGTransport(
    const ParametricMesh& smesh, const KokkosMeshData& kokkosMeshDevice)
    : DGTransport<DG>(smesh)
    , meshDevice(kokkosMeshDevice)
{
    // todo: initialize without hostTransport
    velX = makeKokkosDeviceView("velX", this->velx);
    velY = makeKokkosDeviceView("velY", this->vely);

    normalVelX = makeKokkosDeviceView("normalVelX", this->normalvel_X);
    normalVelY = makeKokkosDeviceView("normalVelY", this->normalvel_Y);

    // parametric map
    std::vector<CG2DGMatrix> cG2DGMatrix;
    cG2DGMatrix.resize(smesh.nelements);
    for (size_t dgi = 0; dgi < smesh.nelements; ++dgi) {
        if (smesh.CoordinateSystem == CARTESIAN) {
            cG2DGMatrix[dgi] = ((ParametricTools::massMatrix<DG>(smesh, dgi).inverse()
                                    * (PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise()
                                        * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, dgi).array()
                                            * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array())
                                              .array())
                                          .matrix())
                * PHI<CGDEGREE, GAUSSPOINTS1D(DG)>.transpose());
        } else {
            cG2DGMatrix[dgi]
                = ((ParametricTools::massMatrix<DG>(smesh, dgi).inverse()
                       * (PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise()
                           * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, dgi).array()
                               * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array()
                               * ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(
                                   smesh, dgi)
                                     .row(1)
                                     .array()
                                     .cos())
                                 .array())
                             .matrix())
                    * PHI<CGDEGREE, GAUSSPOINTS1D(DG)>.transpose());
        }
    }
    cG2DGMatrixDevice = makeKokkosDeviceViewMap("cG2DGMatrix", cG2DGMatrix, true);
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
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> leftEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> rightEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> bottomEdgeOfCell(
    const ConstDeviceViewDG<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> topEdgeOfCell(
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
        return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
    } else if constexpr (DG == 3){
        return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
    } else if constexpr (DG == 6) {
    }
    } else {
        return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
    }
}*/

// dG0 (1 in cell, 1 on edge)

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 1> leftEdgeOfCell(
    const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 1> rightEdgeOfCell(
    const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 1> bottomEdgeOfCell(
    const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid));
}
KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 1> topEdgeOfCell(
    const ConstDeviceViewDG<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid));
}

// dG1 (3 in cell, 2 on edge)

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 2> leftEdgeOfCell(
    const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 2> rightEdgeOfCell(
    const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 2> bottomEdgeOfCell(
    const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 2> topEdgeOfCell(
    const ConstDeviceViewDG<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> leftEdgeOfCell(
    const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5), cv(eid, 4));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> rightEdgeOfCell(
    const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5), cv(eid, 4));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> bottomEdgeOfCell(
    const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5), cv(eid, 3));
}
KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> topEdgeOfCell(
    const ConstDeviceViewDG<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5), cv(eid, 3));
}

// dG2+ (8 in cell, 3 on edge)

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> leftEdgeOfCell(
    const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) - 0.5 * cv(eid, 7));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> rightEdgeOfCell(
    const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) + 0.5 * cv(eid, 7));
}

KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> bottomEdgeOfCell(
    const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) - 0.5 * cv(eid, 6));
}
KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 3> topEdgeOfCell(
    const ConstDeviceViewDG<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) + 0.5 * cv(eid, 6));
}

// add local vector to global row atomically
template <typename T, int Comps, class... Properties>
KOKKOS_IMPL_FUNCTION static void addRowAtomic(const Kokkos::View<T* [Comps], Properties...>& dst,
    const Eigen::Matrix<T, 1, Comps>& src, DeviceIndex row)
{
    for (int i = 0; i < Comps; ++i) {
        Kokkos::atomic_add(&dst(row, i), src(i));
    }
}
// special case because a 1D kokkos view expects only one index
template <typename T, class... Properties>
KOKKOS_IMPL_FUNCTION static void addRowAtomic(
    const Kokkos::View<T*, Properties...>& dst, const Eigen::Matrix<T, 1, 1>& src, DeviceIndex row)
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
            const Eigen::Matrix<FloatType, 1, 2> tangentLeft
                = mesh.edgeVector(ey, ey + mesh.nx + 1);
            const Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> vel1 = 0.5
                * (tangentLeft(0, 1) * leftEdgeOfCell(velXDevice, cy)
                    - tangentLeft(0, 0) * leftEdgeOfCell(velYDevice, cy));
            addRowAtomic(normalVelYDevice, vel1, ey);

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<FloatType, 1, 2> tangentRight
                = mesh.edgeVector(ey + 1, ey + mesh.nx + 2);
            const Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> vel2 = 0.5
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
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangentBottom
                = mesh.edgeVector(nx, nx + 1);
            const Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> vel1 = 0.5
                * (-tangentBottom(0, 1) * bottomEdgeOfCell(velXDevice, cx)
                    + tangentBottom(0, 0) * bottomEdgeOfCell(velYDevice, cx));
            addRowAtomic(normalVelXDevice, vel1, cx);

            /*    normalVelX.row(cx) += 0.5
                    * (-tangentBottom(0, 1) * bottomEdgeOfCell(velXDevice, cx)
                        + tangentBottom(0, 0) * bottomEdgeOfCell(velYDevice, cx));*/

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangentTop
                = mesh.edgeVector(nx + mesh.nx + 1, nx + mesh.nx + 2);
            const Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> vel2 = 0.5
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

template <int DG>
template <int CG>
void KokkosDGTransport<DG>::prepareAdvection(const CGVector<CG>& cgU, const CGVector<CG>& cgV,
    const KokkosDeviceView<CGVector<CG>>& cgUDevice,
    const KokkosDeviceView<CGVector<CG>>& cgVDevice)
{
    DGTransport<DG>::prepareAdvection(cgU, cgV);
/*    auto velXHost = makeKokkosHostView(this->velx);
    auto velYHost = makeKokkosHostView(this->vely);
    Kokkos::deep_copy(velX, velXHost);
    Kokkos::deep_copy(velY, velYHost);*/
    cG2DGDevice(cG2DGMatrixDevice, meshDevice.nx, meshDevice.ny, velX, cgUDevice);
    cG2DGDevice(cG2DGMatrixDevice, meshDevice.nx, meshDevice.ny, velY, cgVDevice);
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

template <int DG>
void KokkosDGTransport<DG>::cG2DGDevice(const KokkosDeviceMapView<CG2DGMatrix>& cG2DGMatrixDevice,
    DeviceIndex nx, DeviceIndex ny, const DeviceViewDG& dgDevice,
    const ConstKokkosDeviceView<CGVector<CGDEGREE>>& cgDevice)
{
    constexpr int CG = CGDEGREE;
//    assert((CG * nx + 1) * (CG * ny + 1) == cg.rows());
//    assert(nx * ny == dg.rows());

    const int cgshift = CG * nx + 1; //!< Index shift for each row

    Kokkos::parallel_for(
        "CG2DG", nx * ny, KOKKOS_LAMBDA(const DeviceIndex dgi) {
            const DeviceIndex iy = dgi / nx; //!< y-index of element
            const DeviceIndex ix = dgi % nx; //!< x-index of element
            const DeviceIndex cgi = CG * cgshift * iy + CG * ix; //!< lower/left Index in cg vector

            Eigen::Matrix<FloatType, (CG == 2 ? 9 : 4), 1>
                cgLocal; //!< the local unknowns in the element
            // we need to use cgDevice outside of the constexpr branch for implicit capture
            const auto& cg = cgDevice;
            if constexpr (CG == 1) {
                cgLocal << cg(cgi), cg(cgi + 1), cg(cgi + cgshift), cg(cgi + 1 + cgshift);
            } else {
                cgLocal << cg(cgi), cg(cgi + 1), cg(cgi + 2), cg(cgi + cgshift),
                    cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift), cg(cgi + 2 * cgshift),
                    cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift);
            }
            // solve:  (Vdg, PHI) = (Vcg, PHI) with mapping to spher. coord.
            auto dg = makeEigenMap(dgDevice);
            dg.row(dgi) = cG2DGMatrixDevice[dgi] * cgLocal;
        });
}

template class KokkosDGTransport<1>;
template class KokkosDGTransport<3>;
template class KokkosDGTransport<6>;
template class KokkosDGTransport<8>;

template void KokkosDGTransport<1>::prepareAdvection(
    const CGVector<CGDEGREE>& cgU, const CGVector<CGDEGREE>& cgV, const KokkosDeviceView<CGVector<CGDEGREE>>& cgUDevice, const KokkosDeviceView<CGVector<CGDEGREE>>& cgVDevice);
template void KokkosDGTransport<3>::prepareAdvection(
    const CGVector<CGDEGREE>& cgU, const CGVector<CGDEGREE>& cgV, const KokkosDeviceView<CGVector<CGDEGREE>>& cgUDevice, const KokkosDeviceView<CGVector<CGDEGREE>>& cgVDevice);
template void KokkosDGTransport<6>::prepareAdvection(
    const CGVector<CGDEGREE>& cgU, const CGVector<CGDEGREE>& cgV, const KokkosDeviceView<CGVector<CGDEGREE>>& cgUDevice, const KokkosDeviceView<CGVector<CGDEGREE>>& cgVDevice);
template void KokkosDGTransport<8>::prepareAdvection(
    const CGVector<CGDEGREE>& cgU, const CGVector<CGDEGREE>& cgV, const KokkosDeviceView<CGVector<CGDEGREE>>& cgUDevice, const KokkosDeviceView<CGVector<CGDEGREE>>& cgVDevice);
}