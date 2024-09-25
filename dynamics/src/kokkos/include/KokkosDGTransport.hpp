/*!
 * @file    KokkosDGTransport.hpp
 * @date    September 12 2024
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSDGTRANSPORT_HPP
#define __KOKKOSDGTRANSPORT_HPP

#include "../include/DGTransport.hpp"
#include "KokkosMeshData.hpp"

namespace Nextsim {
enum struct TimeSteppingScheme { RK1, RK2, RK3 };

template <int DG> constexpr int EDGE_DOFS = (DG == 1) ? 1 : ((DG == 3) ? 2 : 3);

template <int DG> class KokkosDGTransport : DGTransport<DG> {
public:
    using DeviceViewDG = KokkosDeviceView<DGVector<DG>>;
    using HostViewDG = KokkosHostView<DGVector<DG>>;
    using ConstDeviceViewDG = ConstKokkosDeviceView<DGVector<DG>>;

    using DeviceViewEdge = KokkosDeviceView<EdgeVector<EDGE_DOFS<DG>>>;
    using ConstDeviceViewEdge = ConstKokkosDeviceView<EdgeVector<EDGE_DOFS<DG>>>;

    // parametric map
    using CG2DGMatrix = Eigen::Matrix<FloatType, DG, CGDEGREE == 2 ? 9 : 4>;

    KokkosDGTransport(const ParametricMesh& smesh, const KokkosMeshData& kokkosMeshDevice);
    /*!
     * Sets the normal-velocity vector on the edges.
     * The normal velocity is scaled with the length of the edge.
     * This already serves as the integration weight.
     */
    void reinitNormalVelocity();

    /*!
     * Prepares the advection step:
     * - interpolates CG velocity to DG
     * - initializes normal velocity on the edges
     */
    template <int CG> void prepareAdvection(const CGVector<CG>& cgU, const CGVector<CG>& cgV,
        const KokkosDeviceView<CGVector<CG>>& cgUDevice, const KokkosDeviceView<CGVector<CG>>& cgVDevice);
    void step(FloatType dt);

    static void cG2DGDevice(const KokkosDeviceMapView<CG2DGMatrix>& cG2DGMatrixDevice,
        DeviceIndex nx, DeviceIndex ny, const DeviceViewDG& dg,
        const ConstKokkosDeviceView<CGVector<CGDEGREE>>& cg);

    static void reinitNormalVelocityDevice(const DeviceViewEdge& normalVelXDevice,
        const DeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const KokkosMeshData& meshDevice);

private:
    const KokkosMeshData& meshDevice;
    TimeSteppingScheme timeSteppingScheme;

    // current velocity
    DeviceViewDG velX;
    DeviceViewDG velY;

    // normal velocity in edges parallel to X- and Y-axis
    DeviceViewEdge normalVelX;
    DeviceViewEdge normalVelY;

    // temporary vectors for time stepping
    DeviceViewDG tmp1;
    DeviceViewDG tmp2;
    DeviceViewDG tmp3;

    // precomputed maps from ParametricTransportMap
    // KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, GAUSSPOINTS(DG)>> advectionCellTermXDevice;
    // KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, GAUSSPOINTS(DG)>> advectionCellTermYDevice;

    //! The inverse of the dG mass matrix
    // KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>> inverseDGMassMatrixDevice;
    KokkosDeviceMapView<CG2DGMatrix> cG2DGMatrixDevice;
};
}

#endif /* __DGTRANSPORT_HPP */