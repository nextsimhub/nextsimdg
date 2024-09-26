/*!
 * @file    KokkosDGTransport.hpp
 * @date    September 12 2024
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSDGTRANSPORT_HPP
#define __KOKKOSDGTRANSPORT_HPP

#include "../include/CGDynamicsKernel.hpp" // for degree defines
#include "../include/DGTransport.hpp"
#include "KokkosInterpolations.hpp"
#include "KokkosMeshData.hpp"

namespace Nextsim {
enum struct TimeSteppingScheme { RK1, RK2, RK3 };

template <int DG> constexpr int EDGE_DOFS_DG = (DG == 1) ? 1 : ((DG == 3) ? 2 : 3);
template <int DG>
constexpr int GAUSS_POINTS = ((DG == 8) || (DG == 6)) ? 9
    : (DG == 3)                                       ? 4
    : (DG == 1)                                       ? 1
                                                      : -1;

template <int DG> class KokkosDGTransport : DGTransport<DG> {
public:
    using DeviceViewDG = KokkosDeviceView<DGVector<DG>>;
    using HostViewDG = KokkosHostView<DGVector<DG>>;
    using ConstDeviceViewDG = ConstKokkosDeviceView<DGVector<DG>>;

    static constexpr int EDGE_DOFS = EDGE_DOFS_DG<DG>;
    using DeviceViewEdge = KokkosDeviceView<EdgeVector<EDGE_DOFS>>;
    using ConstDeviceViewEdge = ConstKokkosDeviceView<EdgeVector<EDGE_DOFS>>;

    static constexpr int GP = GAUSS_POINTS<DG>;
    static constexpr int GP1D = GAUSSPOINTS1D(DG);
    using AdvectionCellTerm = Eigen::Matrix<FloatType, DG, GP>;

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

    void step(FloatType dt, DGVector<DG>& phi);

    static void reinitNormalVelocityDevice(const DeviceViewEdge& normalVelXDevice,
        const DeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const KokkosMeshData& meshDevice);

    void addCellTermsDevice(FloatType dt, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const ConstDeviceViewDG& phiDevice,
        const DeviceViewDG& phiupDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
        const ConstDeviceBitset& landMaskDevice);

    void addEdgeXTermsDevice(FloatType dt, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& phiDevice,
        const DeviceViewDG& phiupDevice, const ConstDeviceBitset& landMaskDevice, DeviceIndex nx,
        DeviceIndex ny);

    void dGTransportOperatorDevice(FloatType dt, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& phiDevice,
        const DeviceViewDG& phiupDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
        const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
        const KokkosMeshData& meshDevice);

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
    DeviceViewDG tmpRes1;
    DeviceViewDG tmpRes2;
    DeviceViewDG tmpRes3;

    // precomputed maps from ParametricTransportMap
    KokkosDeviceMapView<AdvectionCellTerm> advectionCellTermXDevice;
    KokkosDeviceMapView<AdvectionCellTerm> advectionCellTermYDevice;

    //! The inverse of the dG mass matrix
    KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>> inverseDGMassMatrixDevice;

    // using PSI1DView = ConstKokkosDeviceView<decltype(PSI<DG, GAUSSPOINTS1D(DG)>)>;
    // using PSIeView = ConstKokkosDeviceView<decltype(PSI<DG, GAUSSPOINTS1D(DG)>)>;
    // using PSIewView = ConstKokkosDeviceView<decltype(PSI<DG, GAUSSPOINTS1D(DG)>)>;
};
}

#endif /* __DGTRANSPORT_HPP */