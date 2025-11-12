/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

// by also guarding for USE_KOKKOS this header can be safely included even
// when Kokkos is not enabled
#if !defined(KOKKOSMEVPDYNAMICSKERNEL_HPP) && defined(USE_KOKKOS)
#define KOKKOSMEVPDYNAMICSKERNEL_HPP

#include "../../include/VPParameters.hpp"
#include "KokkosCGDynamicsKernel.hpp"
#include "../../../core/src/kokkos/include/KokkosUtils.hpp"

namespace Nextsim {

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection>
class KokkosMEVPDynamicsKernel : public KokkosCGDynamicsKernel<DGadvection> {
private:
    using EdgeVec = typename KokkosCGDynamicsKernel<DGadvection>::EdgeVec;

public:
    using Base = KokkosCGDynamicsKernel<DGadvection>;

    using DeviceViewCG = typename Base::DeviceViewCG;
    using ConstDeviceViewCG = typename Base::ConstDeviceViewCG;

    using DeviceViewAdvect = typename Base::DeviceViewAdvect;
    using HostViewAdvect = typename Base::HostViewAdvect;
    using ConstDeviceViewAdvect = typename Base::ConstDeviceViewAdvect;

    using DeviceViewStress = typename Base::DeviceViewStress;
    using ConstDeviceViewStress = typename Base::ConstDeviceViewStress;

    using PSIAdvectView = typename Base::PSIAdvectView;
    using PSIStressView = typename Base::PSIStressView;

    using GaussMapDevice = typename Base::GaussMapDevice;

    KokkosMEVPDynamicsKernel(const VPParameters& paramsIn);

    KokkosMEVPDynamicsKernel(const KokkosMEVPDynamicsKernel<DGadvection>&) = delete;
    KokkosMEVPDynamicsKernel(KokkosMEVPDynamicsKernel<DGadvection>&&) = delete;

    KokkosMEVPDynamicsKernel& operator=(const KokkosMEVPDynamicsKernel<DGadvection>&) = delete;
    KokkosMEVPDynamicsKernel& operator=(KokkosMEVPDynamicsKernel<DGadvection>&&) = delete;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void update(const TimestepTime& tst) override;

    // cuda requires these functions to be public
    static void updateStressHighOrderDevice(const DeviceViewStress& s11Device,
        const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const PSIAdvectView& PSIAdvectDevice,
        const PSIStressView& PSIStressDevice, const ConstDeviceViewAdvect& hiceDevice,
        const ConstDeviceViewAdvect& ciceDevice, const ConstDeviceBitset& landMaskDevice,
        const GaussMapDevice& iMJwPSIDevice, const VPParameters& params, FloatType alpha);
    static void updateMomentumDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const ConstDeviceViewCG& u0Device, const ConstDeviceViewCG& v0Device,
        const ConstDeviceViewCG& cgHDevice, const ConstDeviceViewCG& cgADevice,
        const ConstDeviceViewCG& uAtmosDevice, const ConstDeviceViewCG& vAtmosDevice,
        const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
        const ConstDeviceViewCG& dStressXDevice, const ConstDeviceViewCG& dStressYDevice,
        const ConstDeviceViewCG& xGradSeaSurfaceHeightDevice,
        const ConstDeviceViewCG& yGradSeaSurfaceHeightDevice,
        const ConstDeviceViewCG& lumpedCGMassDevice, const ConstDeviceBitset& cgLandMaskDevice,
        const TimestepTime& tst, const VPParameters& params, FloatType beta);

private:
    // Step-initial ice velocity
    // mutable variants are needed to copy the data but accesses on the device are all constant
    DeviceViewCG u0DeviceMut;
    DeviceViewCG v0DeviceMut;
    ConstDeviceViewCG u0Device;
    ConstDeviceViewCG v0Device;

    const VPParameters& params;
    FloatType alpha = 1500.;
    FloatType beta = 1500.;
};

} /* namespace Nextsim */

#endif /* KOKKOSMEVPDYNAMICSKERNEL_HPP */
