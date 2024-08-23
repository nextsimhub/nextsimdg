/*!
 * @file KokkosCGDynamicsKernel.hpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSCGDYNAMICSKERNEL_HPP
#define KOKKOSCGDYNAMICSKERNEL_HPP

#include "../../include/CGDynamicsKernel.hpp"
#include "KokkosMeshData.hpp"

namespace Nextsim {

template <int DG> constexpr int NGP_DG = ((DG == 8) || (DG == 6)) ? 3 : (DG == 3 ? 2 : -1);

using DeviceIndex = EIGEN_DEFAULT_DENSE_INDEX_TYPE;

template <int DGadvection> class KokkosCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
public:
    // common types for Kokkos buffers
    // cG components
    using DeviceViewCG = KokkosDeviceView<CGVector<CGdegree>>;
    using HostViewCG = KokkosHostView<CGVector<CGdegree>>;
    using ConstDeviceViewCG = ConstKokkosDeviceView<CGVector<CGdegree>>;

    // strain and stress components
    using DeviceViewStress = KokkosDeviceView<DGVector<DGstressComp>>;
    using HostViewStress = KokkosHostView<DGVector<DGstressComp>>;
    using ConstDeviceViewStress = ConstKokkosDeviceView<DGVector<DGstressComp>>;
    // using DeviceSymmetricTensorVector = std::array<DeviceViewStress, 3>;

    // precomputed maps
    using DivMapDevice = KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix>;
    using GradMapDevice = KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix>;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    static void dirichletZero(const DeviceViewCG& v, DeviceIndex nx, DeviceIndex ny,
        const KokkosMeshData::DirichletData& dirichlet);
    // cuda requires these functions to be public but they should only be needed by the concrete
    // dynamics (like protected)
    static void projectVelocityToStrainDevice(const ConstDeviceViewCG& uDevice,
        const ConstDeviceViewCG& vDevice, const DeviceViewStress& e11Device,
        const DeviceViewStress& e12Device, const DeviceViewStress& e22Device,
        const ConstDeviceBitset& landMaskDevice, const GradMapDevice& iMgradXDevice,
        const GradMapDevice& iMgradYDevice, const GradMapDevice& iMMDevice, DeviceIndex nx,
        DeviceIndex ny, COORDINATES coordinates);
    static void computeStressDivergenceDevice(const DeviceViewCG& dStressXDevice,
        const DeviceViewCG& dStressYDevice, const ConstDeviceViewStress& s11Device,
        const ConstDeviceViewStress& s12Device, const ConstDeviceViewStress& s22Device,
        const ConstDeviceBitset& landMaskDevice, const DivMapDevice& divS1Device,
        const DivMapDevice& divS2Device, const DivMapDevice& divMDevice,
        const KokkosMeshData::DirichletData& dirichletDevice, DeviceIndex nx, DeviceIndex ny,
        COORDINATES coordinates);
    static void applyBoundariesDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const KokkosMeshData::DirichletData& dirichletDevice, DeviceIndex nx, DeviceIndex ny);

protected:
    // cG (velocity) components
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

    // dG stress
    DeviceViewStress s11Device;
    HostViewStress s11Host;
    DeviceViewStress s12Device;
    HostViewStress s12Host;
    DeviceViewStress s22Device;
    HostViewStress s22Host;
    DeviceViewStress e11Device;
    HostViewStress e11Host;
    DeviceViewStress e12Device;
    HostViewStress e12Host;
    DeviceViewStress e22Device;
    HostViewStress e22Host;

    // precomputed parametric map
    DivMapDevice divS1Device;
    DivMapDevice divS2Device;
    DivMapDevice divMDevice;

    GradMapDevice iMgradXDevice;
    GradMapDevice iMgradYDevice;
    GradMapDevice iMMDevice;

    std::unique_ptr<KokkosMeshData> meshData;
};

}

#endif // KOKKOSCGDYNAMICSKERNEL_HPP