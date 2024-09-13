/*!
 * @file KokkosDGTransport.cpp
 * @date September 12, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosDGTransport.hpp"

namespace Nextsim {

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
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> leftEdgeOfCell(const DGVector<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> rightEdgeOfCell(const DGVector<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> bottomEdgeOfCell(
    const DGVector<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGE_DOFS<DG>> topEdgeOfCell(const DGVector<DG>& cv, DeviceIndex eid);

namespace Details {
    template <typename DGVec> struct DGDegree;

    template <typename T, int DG, class... Properties>
    struct DGDegree<Kokkos::View<T* [DG], Properties...>> {
        static constexpr int value = DG;
    };
}

template <typename DGVec>
static KOKKOS_IMPL_FUNCTION auto leftEdgeOfCell(const DGVec& cv, DeviceIndex eid)
{
    constexpr int DG = Details::DGDegree<DGVec>::value;
    if constexpr (DG == 1) {
        return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
    } else {
        return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
    }
}

// dG0 (1 in cell, 1 on edge)
/*
Eigen::Matrix<FloatType, 1, 1> leftEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}*/

Eigen::Matrix<FloatType, 1, 1> rightEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}

Eigen::Matrix<FloatType, 1, 1> bottomEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}
Eigen::Matrix<FloatType, 1, 1> topEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}

// dG1 (3 in cell, 2 on edge)

Eigen::Matrix<FloatType, 1, 2> leftEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}

Eigen::Matrix<FloatType, 1, 2> rightEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}

Eigen::Matrix<FloatType, 1, 2> bottomEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
Eigen::Matrix<FloatType, 1, 2> topEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)

Eigen::Matrix<FloatType, 1, 3> leftEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5), cv(eid, 4));
}

Eigen::Matrix<FloatType, 1, 3> rightEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5), cv(eid, 4));
}

Eigen::Matrix<FloatType, 1, 3> bottomEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5), cv(eid, 3));
}
Eigen::Matrix<FloatType, 1, 3> topEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5), cv(eid, 3));
}

// dG2+ (8 in cell, 3 on edge)

Eigen::Matrix<FloatType, 1, 3> leftEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) - 0.5 * cv(eid, 7));
}

Eigen::Matrix<FloatType, 1, 3> rightEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) + 0.5 * cv(eid, 7));
}

Eigen::Matrix<FloatType, 1, 3> bottomEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) - 0.5 * cv(eid, 6));
}
Eigen::Matrix<FloatType, 1, 3> topEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) + 0.5 * cv(eid, 6));
}

template <int DG>
void KokkosDGTransport<DG>::reinitNormalVelocityDevice(const DeviceViewEdge& normalVelX,
    const DeviceViewEdge& normalVelY, const ConstDeviceViewDG& velX, const ConstDeviceViewDG& velY,
    const ConstDeviceBitset& landMaskDevice, DeviceIndex nx, DeviceIndex ny,
    const KokkosMeshData::DirichletData& dirichletDevice)
{
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, normalVelX, 0.0);
    Kokkos::deep_copy(execSpace, normalVelY, 0.0);

 /*   // average the velocity to the Y-edges
#pragma omp parallel for
    for (DeviceIndex iy = 0; iy < ny; ++iy) {
        //   |     |
        // --*-----*--
        //  ey  cy |
        //   |     |
        // -ey-----*--
        //   |     |
        DeviceIndex ey = iy * (nx + 1); // first edge-index and node-index in row
        DeviceIndex cy = iy * nx; // first cell index in row

        for (DeviceIndex ix = 0; ix < nx; ++ix, ++ey, ++cy) {
            if (!landMaskDevice.test(cy)) {
                return;
            }

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangentLeft
                = smesh.edgevector(ey, ey + nx + 1);
            normalvel_Y.row(ey) += 0.5
                * (tangentLeft(0, 1) * leftedgeofcell<DG>(velx, cy)
                    - tangentLeft(0, 0) * leftedgeofcell<DG>(vely, cy));

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangentRight
                = smesh.edgevector(ey + 1, ey + nx + 2);
            normalvel_Y.row(ey + 1) += 0.5
                * (tangentRight(0, 1) * rightedgeofcell<DG>(velx, cy)
                    - tangentRight(0, 0) * rightedgeofcell<DG>(vely, cy));
        }
        // we need an adjustment along the boundaries.. This is done later on.
    }

#pragma omp parallel for
    for (DeviceIndex ix = 0; ix < nx; ++ix) {
        //   |     |
        // --*-----*--
        //   |  cx |
        //   |     |
        // -nx-ex--*--
        //   |     |

        DeviceIndex cx = ix; // first edge-index and cell-index
        DeviceIndex nx = ix; // first cell index in row

        for (DeviceIndex iy = 0; iy < ny; ++iy, cx += nx, nx += nx + 1) {
            if (smesh.landmask[cx] == 0) // skip land elements
                continue;

            // un-normed tangent vector of bottom edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_bottom
                = smesh.edgevector(nx, nx + 1);

            normalvel_X.row(cx) += 0.5
                * (-tangent_bottom(0, 1) * bottomedgeofcell<DG>(velx, cx)
                    + tangent_bottom(0, 0) * bottomedgeofcell<DG>(vely, cx));

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_top
                = smesh.edgevector(nx + nx + 1, nx + nx + 2);

            normalvel_X.row(cx + nx) += 0.5
                * (-tangent_top(0, 1) * topedgeofcell<DG>(velx, cx)
                    + tangent_top(0, 0) * topedgeofcell<DG>(vely, cx));
        }
    }

    // Take care of the boundaries. Usually, the normal velocity is the average velocity
    // from left and from the right. Hence, we get the factor 0.5 above. At boundaries,
    // the normal is set only once, from the inside. These edges must be scaled with 2.0

    for (DeviceIndex seg = 0; seg < 4; ++seg) // run over the 4 segments (bot, right, top, left)
    {
#pragma omp parallel for
        for (DeviceIndex i = 0; i < smesh.dirichlet[seg].size(); ++i) {
            const DeviceIndex eid = smesh.dirichlet[seg][i]; //! The i of the boundary element
            const DeviceIndex ix = eid % nx; //! x & y indices of the element
            const DeviceIndex iy = eid / nx;

            if (seg == 0) // bottom
                normalvel_X.row(nx * iy + ix) *= 2.0;
            else if (seg == 1) // right
                normalvel_Y.row((nx + 1) * iy + ix + 1) *= 2.0;
            else if (seg == 2) // top
                normalvel_X.row(nx * (iy + 1) + ix) *= 2.0;
            else if (seg == 3) // left
                normalvel_Y.row((nx + 1) * iy + ix) *= 2.0;
        }
    }*/
}

}