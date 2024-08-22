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

#include "../../include/CGDynamicsKernel.hpp"
#include "../../include/VPParameters.hpp"
#include "KokkosUtils.hpp"

#include <Kokkos_Bitset.hpp>

namespace Nextsim {

template <int DG> constexpr int NGP_DG = ((DG == 8) || (DG == 6)) ? 3 : (DG == 3 ? 2 : -1);
// todo: move into kokkos / gpu namespace
using DeviceIndex = EIGEN_DEFAULT_DENSE_INDEX_TYPE;

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection> class KokkosMEVPDynamicsKernel : public CGDynamicsKernel<DGadvection> {
private:
    static constexpr int NGP = NGP_DG<DGstressComp>;

    using EdgeVec = Eigen::Matrix<FloatType, 1, NGP * NGP>;

public:
    struct KokkosBuffers {
        // cG (velocity) components
        using DeviceViewCG = KokkosDeviceView<CGVector<CGdegree>>;
        using HostViewCG = KokkosHostView<CGVector<CGdegree>>;
        using ConstDeviceViewCG = ConstKokkosDeviceView<CGVector<CGdegree>>;
        DeviceViewCG uDevice;
        HostViewCG uHost;
        DeviceViewCG vDevice;
        HostViewCG vHost;
        // Step-initial ice velocity
        // mutable variants are needed to copy the data but accesses on the device are all constant
        DeviceViewCG u0DeviceMut;
        DeviceViewCG v0DeviceMut;
        ConstDeviceViewCG u0Device;
        ConstDeviceViewCG v0Device;

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

        // strain and stress components
        using DeviceViewStress = KokkosDeviceView<DGVector<DGstressComp>>;
        using HostViewStress = KokkosHostView<DGVector<DGstressComp>>;
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

        using DeviceViewAdvect = KokkosDeviceView<DGVector<DGadvection>>;
        using HostViewAdvect = KokkosHostView<DGVector<DGadvection>>;
        DeviceViewAdvect hiceDevice;
        HostViewAdvect hiceHost;
        DeviceViewAdvect ciceDevice;
        HostViewAdvect ciceHost;

        DeviceViewCG cgHDevice;
        HostViewCG cgHHost;
        DeviceViewCG cgADevice;
        HostViewCG cgAHost;

        // constant matrices also need to be available on the GPU
        using PSIAdvectType = decltype(PSI<DGadvection, NGP>);
        using PSIStressType = decltype(PSI<DGstressComp, NGP>);
        ConstKokkosDeviceView<PSIAdvectType> PSIAdvectDevice;
        ConstKokkosDeviceView<PSIStressType> PSIStressDevice;

        // parametric map precomputed transforms
        // todo: refactor into KokkosParametricMap with switch for precomputed / on-the-fly
        ConstDeviceViewCG lumpedcgmassDevice;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix> divS1Device;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix> divS2Device;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix> divMDevice;

        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix> iMgradXDevice;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix> iMgradYDevice;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix> iMMDevice;

        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix> iMJwPSIDevice;

        // mesh related
        std::array<KokkosDeviceMapView<size_t>, 4> dirichletDevice;
        Kokkos::ConstBitset<Kokkos::DefaultExecutionSpace> landMaskDevice;
    };

    KokkosMEVPDynamicsKernel(const VPParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>()
        , params(paramsIn)
    {
    }

    KokkosMEVPDynamicsKernel(const KokkosMEVPDynamicsKernel<DGadvection>&) = delete;
    KokkosMEVPDynamicsKernel(KokkosMEVPDynamicsKernel<DGadvection>&&) = delete;

    KokkosMEVPDynamicsKernel& operator=(const KokkosMEVPDynamicsKernel<DGadvection>&) = delete;
    KokkosMEVPDynamicsKernel& operator=(KokkosMEVPDynamicsKernel<DGadvection>&&) = delete;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void update(const TimestepTime& tst) override;

    // todo: move kokkos stuff out of class into extra namespace?
    // todo: rename functions that currently conflict with base class versions
    // cuda requires these functions to be public
    static void projVelocityToStrain(
        const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates);
    static void updateStressHighOrder(
        const KokkosBuffers& _buffers, const VPParameters& _params, FloatType _alpha);
    static void computeStressDivergence(
        const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates);
    static void applyBoundariesDevice(
        const KokkosBuffers& _buffers, DeviceIndex nx, DeviceIndex ny);
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
