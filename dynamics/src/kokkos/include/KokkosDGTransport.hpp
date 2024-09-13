/*!
 * @file    KokkosDGTransport.hpp
 * @date    September 12 2024
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSDGTRANSPORT_HPP
#define __KOKKOSDGTRANSPORT_HPP

#include "KokkosMeshData.hpp"
#include "../include/DGTransport.hpp"

namespace Nextsim {
enum struct TimeSteppingScheme { RK1, RK2, RK3 };

template<int DG>
constexpr int EDGE_DOFS = (DG == 1) ? 1 : ((DG == 3) ? 2 : 3);

template <int DG> class KokkosDGTransport {
public:
    using DeviceViewDG = KokkosDeviceView<DGVector<DG>>;
    using HostViewDG = KokkosHostView<DGVector<DG>>;
    using ConstDeviceViewDG = ConstKokkosDeviceView<DGVector<DG>>;

    using DeviceViewEdge = KokkosDeviceView<EdgeVector<EDGE_DOFS<DG>>>;
    using ConstDeviceViewEdge = ConstKokkosDeviceView<EdgeVector<EDGE_DOFS<DG>>>;

    /*!
     * Sets the normal-velocity vector on the edges
     * The normal velocity is scaled with the length of the edge,
     * this already serves as the integraiton weight
     */
    void reinitNormalVelocity();

    /*!
     * Prepares the advection step:
     * - interpolates CG velocity to DG
     * - initializes normal velocity on the edges
     */
    template <int CG> void prepareAdvection(const CGVector<CG>& cgU, const CGVector<CG>& cgV);
    void step(FloatType dt);

    static void reinitNormalVelocityDevice(const DeviceViewEdge& normalVelX,
        const DeviceViewEdge& normalVelY, const ConstDeviceViewDG& velX,
        const ConstDeviceViewDG& velY, const ConstDeviceBitset& landMaskDevice, DeviceIndex nx,
        DeviceIndex ny, const KokkosMeshData::DirichletData& dirichletDevice);

private:
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
};
}

#endif /* __DGTRANSPORT_HPP */