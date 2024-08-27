/*!
 * @file KokkosMEVPDynamicsKernel.hpp
 *
 * @date Feb 2, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

// by also guarding for USE_KOKKOS this header can be safely included even
// when Kokkos is not enabled
#if not defined KOKKOSMEVPDYNAMICSKERNEL_HPP && USE_KOKKOS
#define KOKKOSMEVPDYNAMICSKERNEL_HPP

#include "../../include/VPParameters.hpp"
#include "KokkosCGDynamicsKernel.hpp"
#include "KokkosUtils.hpp"

namespace Nextsim {

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection>
class KokkosMEVPDynamicsKernel : public KokkosCGDynamicsKernel<DGadvection> {
private:
    static constexpr int NGP = NGP_DG<DGstressComp>;

    using EdgeVec = Eigen::Matrix<FloatType, 1, NGP * NGP>;

public:
    using Base = KokkosCGDynamicsKernel<DGadvection>;
    using DeviceViewCG = typename Base::DeviceViewCG;
    using ConstDeviceViewCG = typename Base::ConstDeviceViewCG;

    using DeviceViewAdvect = KokkosDeviceView<DGVector<DGadvection>>;
    using HostViewAdvect = KokkosHostView<DGVector<DGadvection>>;
    using ConstDeviceViewAdvect = ConstKokkosDeviceView<DGVector<DGadvection>>;

    using DeviceViewStress = typename Base::DeviceViewStress;
    using ConstDeviceViewStress = typename Base::ConstDeviceViewStress;

    using PSIAdvectType = decltype(PSI<DGadvection, NGP>);
    using PSIStressType = decltype(PSI<DGstressComp, NGP>);
    // in gcc13 the signature of updateStressHighOrder is somehow incompatible between declaration
    // and implementation if we use ConstKokkosDeviceView directly
    using PSIAdvectView = ConstKokkosDeviceView<PSIAdvectType>;
    using PSIStressView = ConstKokkosDeviceView<PSIStressType>;

    KokkosMEVPDynamicsKernel(const VPParameters& paramsIn)
        : KokkosCGDynamicsKernel<DGadvection>()
        , params(paramsIn)
    {
    }

    KokkosMEVPDynamicsKernel(const KokkosMEVPDynamicsKernel<DGadvection>&) = delete;
    KokkosMEVPDynamicsKernel(KokkosMEVPDynamicsKernel<DGadvection>&&) = delete;

    KokkosMEVPDynamicsKernel& operator=(const KokkosMEVPDynamicsKernel<DGadvection>&) = delete;
    KokkosMEVPDynamicsKernel& operator=(KokkosMEVPDynamicsKernel<DGadvection>&&) = delete;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void update(const TimestepTime& tst) override;

    // cuda requires these functions to be public
    static void updateStressHighOrder(const DeviceViewStress& s11Device,
        const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const PSIAdvectView& PSIAdvectDevice,
        const PSIStressView& PSIStressDevice, const ConstDeviceViewAdvect& hiceDevice,
        const ConstDeviceViewAdvect& ciceDevice,
        const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix>& iMJwPSIDevice,
        const VPParameters& params, FloatType alphaa);
    static void updateMomentumDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const ConstDeviceViewCG& u0Device, const ConstDeviceViewCG& v0Device,
        const ConstDeviceViewCG& cgHDevice, const ConstDeviceViewCG& cgADevice,
        const ConstDeviceViewCG& uAtmosDevice, const ConstDeviceViewCG& vAtmosDevice,
        const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
        const ConstDeviceViewCG& dStressXDevice, const ConstDeviceViewCG& dStressYDevice,
        const ConstDeviceViewCG& lumpedCGMassDevice, const TimestepTime& tst,
        const VPParameters& params, FloatType beta);

private:
    // Step-initial ice velocity
    // mutable variants are needed to copy the data but accesses on the device are all constant
    DeviceViewCG u0DeviceMut;
    DeviceViewCG v0DeviceMut;
    ConstDeviceViewCG u0Device;
    ConstDeviceViewCG v0Device;

    DeviceViewAdvect hiceDevice;
    HostViewAdvect hiceHost;
    DeviceViewAdvect ciceDevice;
    HostViewAdvect ciceHost;

    // constant matrices also need to be available on the GPU
    ConstKokkosDeviceView<PSIAdvectType> PSIAdvectDevice;
    ConstKokkosDeviceView<PSIStressType> PSIStressDevice;

    // parametric map precomputed transforms
    // todo: refactor into KokkosParametricMap with switch for precomputed / on-the-fly
    ConstDeviceViewCG lumpedcgmassDevice;
    KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix> iMJwPSIDevice;

    const VPParameters& params;
    FloatType alpha = 1500.;
    FloatType beta = 1500.;

    // todo: change base class to remove this completely
    void updateMomentum(const TimestepTime& tst) override { }
};

} /* namespace Nextsim */

#endif /* KOKKOSMEVPDYNAMICSKERNEL_HPP */
