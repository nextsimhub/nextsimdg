/*!
 * @file KokkosCGDynamicsKernel.hpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSCGDYNAMICSKERNEL_HPP
#define KOKKOSCGDYNAMICSKERNEL_HPP

#include "../../include/CGDynamicsKernel.hpp"

namespace Nextsim {

template <int DG> constexpr int NGP_DG = ((DG == 8) || (DG == 6)) ? 3 : (DG == 3 ? 2 : -1);
// todo: move into kokkos / gpu namespace
using DeviceIndex = EIGEN_DEFAULT_DENSE_INDEX_TYPE;

template <int DGadvection> class KokkosCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
public:
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    // cuda requires these functions to be public but they should only be needed by the concrete
    // dynamics (like protected)
    static void projectVelocityToStrainDevice(
        const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates);
    static void computeStressDivergenceDevice(const KokkosBuffers& _buffers,
        const KokkosMeshData& _meshData, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates);
    static void applyBoundariesDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const std::array<KokkosDeviceMapView, 4>& dirichlet, DeviceIndex nx, DeviceIndex ny);

private:
    // cG (velocity) components
    using DeviceViewCG = KokkosDeviceView<CGVector<CGdegree>>;
    using HostViewCG = KokkosHostView<CGVector<CGdegree>>;
    using ConstDeviceViewCG = ConstKokkosDeviceView<CGVector<CGdegree>>;
    DeviceViewCG uDevice;
    HostViewCG uHost;
    DeviceViewCG vDevice;
    HostViewCG vHost;

    DeviceViewCG cgHDevice;
    HostViewCG cgHHost;
    DeviceViewCG cgADevice;
    HostViewCG cgAHost;

    DeviceViewCG dStressXDevice;
    HostViewCG dStressXHost;
    DeviceViewCG dStressYDevice;
    HostViewCG dStressYHost;

    DeviceViewCG uOceanDevice;
    HostViewCG uOceanHost;
    DeviceViewCG vOceanDevice;
    HostViewCG vOceanHost;

    DeviceViewCG uAtmosDevice;
    HostViewCG uAtmosHost;
    DeviceViewCG vAtmosDevice;
    HostViewCG vAtmosHost;

    std::unique_ptr<KokkosMeshData> meshData;
};

}

#endif // KOKKOSCGDYNAMICSKERNEL_HPP