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

    void setTimeSteppingScheme(TimeSteppingScheme tss);
    /*!
     * Prepares the advection step:
     * - interpolates CG velocity to DG
     * - initializes normal velocity on the edges
     */
    void prepareAdvection(const KokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
        const KokkosDeviceView<CGVector<CGdegree>>& cgVDevice);

    void step(FloatType dt, const DeviceViewDG& phiDevice);

    // internal methods, only public because cuda needs it
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
        const ConstDeviceViewDG& phiDevice, const DeviceViewDG& phiupDevice,
        const ConstDeviceBitset& landMaskDevice, DeviceIndex nx, DeviceIndex ny);

    void addEdgeYTermsDevice(FloatType dt, const ConstDeviceViewEdge& normalVelYDevice,
        const ConstDeviceViewDG& phiDevice, const DeviceViewDG& phiupDevice,
        const ConstDeviceBitset& landMaskDevice, DeviceIndex nx, DeviceIndex ny);

    void addBoundaryTermsDevice(FloatType dt, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& phiDevice,
        const DeviceViewDG& phiupDevice, const KokkosMeshData& meshDevice);

    void dGTransportOperatorDevice(FloatType dt, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const ConstDeviceViewDG& phiDevice,
        const DeviceViewDG& phiupDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
        const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
        const KokkosMeshData& meshDevice);

    void stepRK1(FloatType dt, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
        const DeviceViewDG& tmpRes1,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
        const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
        const KokkosMeshData& meshDevice);

    void stepRK2(FloatType dt, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
        const DeviceViewDG& tmpRes1, const DeviceViewDG& tmpRes2,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermXDevice,
        const KokkosDeviceMapView<AdvectionCellTerm>& advectionCellTermYDevice,
        const KokkosDeviceMapView<Eigen::Matrix<FloatType, DG, DG>>& inverseDGMassMatrixDevice,
        const KokkosMeshData& meshDevice);

    void stepRK3(FloatType dt, const ConstDeviceViewDG& velXDevice,
        const ConstDeviceViewDG& velYDevice, const ConstDeviceViewEdge& normalVelXDevice,
        const ConstDeviceViewEdge& normalVelYDevice, const DeviceViewDG& phiDevice,
        const DeviceViewDG& tmpRes1, const DeviceViewDG& tmpRes2, const DeviceViewDG& tmpRes3,
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