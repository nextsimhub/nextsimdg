/*!
 * @file    KokkosDGTransport.hpp
 * @date    September 12 2024
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSDGTRANSPORT_HPP
#define __KOKKOSDGTRANSPORT_HPP

#include "../include/DGTransport.hpp"
#include "../include/CGDynamicsKernel.hpp" // for degree defines
#include "KokkosMeshData.hpp"
#include "KokkosInterpolations.hpp"

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

    KokkosDGTransport(const ParametricMesh& smesh, const KokkosMeshData& _meshDevice,
        const Interpolations::KokkosCG2DGInterpolator<DG, CGdegree>& _cG2DGInterpolator);
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
    void prepareAdvection(const CGVector<CGdegree>& cgU, const CGVector<CGdegree>& cgV,
        const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
        const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice);
    void step(FloatType dt);

    static void reinitNormalVelocityDevice(const DeviceViewEdge& normalVelXDevice,
        const DeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const KokkosMeshData& meshDevice);

private:
    const KokkosMeshData& meshDevice;
    const Interpolations::KokkosCG2DGInterpolator<DG, CGdegree>& cG2DGInterpolator;
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
};
}

#endif /* __DGTRANSPORT_HPP */