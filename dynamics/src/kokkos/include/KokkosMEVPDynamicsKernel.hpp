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
    struct KokkosBuffers {
        // Step-initial ice velocity
        // mutable variants are needed to copy the data but accesses on the device are all constant
        using DeviceViewCG = typename KokkosCGDynamicsKernel<DGadvection>::DeviceViewCG;
        using ConstDeviceViewCG = typename KokkosCGDynamicsKernel<DGadvection>::ConstDeviceViewCG;
        DeviceViewCG u0DeviceMut;
        DeviceViewCG v0DeviceMut;
        ConstDeviceViewCG u0Device;
        ConstDeviceViewCG v0Device;

        using DeviceViewAdvect = KokkosDeviceView<DGVector<DGadvection>>;
        using HostViewAdvect = KokkosHostView<DGVector<DGadvection>>;
        DeviceViewAdvect hiceDevice;
        HostViewAdvect hiceHost;
        DeviceViewAdvect ciceDevice;
        HostViewAdvect ciceHost;

        // constant matrices also need to be available on the GPU
        using PSIAdvectType = decltype(PSI<DGadvection, NGP>);
        using PSIStressType = decltype(PSI<DGstressComp, NGP>);
        ConstKokkosDeviceView<PSIAdvectType> PSIAdvectDevice;
        ConstKokkosDeviceView<PSIStressType> PSIStressDevice;

        // parametric map precomputed transforms
        // todo: refactor into KokkosParametricMap with switch for precomputed / on-the-fly
        ConstDeviceViewCG lumpedcgmassDevice;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix> iMJwPSIDevice;
    };

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
    static void updateStressHighOrder(
        const KokkosBuffers& _buffers, const VPParameters& _params, FloatType _alpha);
    static void updateMomentumDevice(const TimestepTime& tst, const KokkosBuffers& _buffers,
        const VPParameters& _params, FloatType beta);

private:
    KokkosBuffers buffers;
    const VPParameters& params;
    FloatType alpha = 1500.;
    FloatType beta = 1500.;

    // todo: change base class to remove this completely
    void updateMomentum(const TimestepTime& tst) override { }
};

} /* namespace Nextsim */

#endif /* KOKKOSMEVPDYNAMICSKERNEL_HPP */
