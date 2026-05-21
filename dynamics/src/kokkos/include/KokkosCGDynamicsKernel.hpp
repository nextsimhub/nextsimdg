/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSCGDYNAMICSKERNEL_HPP
#define KOKKOSCGDYNAMICSKERNEL_HPP

#include "../../../core/src/kokkos/include/KokkosModelArray.hpp"
#include "../../../core/src/kokkos/include/KokkosUtils.hpp"
#include "../../include/CGDynamicsKernel.hpp"
#include "KokkosDGModelArray.hpp"

namespace Nextsim {

// predeclarations for pimpl
class KokkosMesh;
template <int DG> class KokkosSlopeLimiter;
template <int DG> class KokkosDGTransport;

namespace Interpolations {
    template <int DG, int CG> class KokkosCG2DGInterpolator;
    template <int CG, int DG> class KokkosDG2CGInterpolator;
}

template <int DG> constexpr int NGP_DG = ((DG == 8) || (DG == 6)) ? 3 : (DG == 3 ? 2 : -1);

template <int DGadvection> class KokkosCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
public:
    using Base = CGDynamicsKernel<DGadvection>;
    // common types for Kokkos buffers
    // cG components
    using DeviceViewCG = KokkosDeviceView<CGVector<CGdegree>>;
    using HostViewCG = KokkosHostView<CGVector<CGdegree>>;
    using ConstDeviceViewCG = ConstKokkosDeviceView<CGVector<CGdegree>>;
    using DeviceViewCG1 = KokkosDeviceView<CGVector<1>>;
    using ConstDeviceViewCG1 = ConstKokkosDeviceView<CGVector<1>>;

    using DeviceViewDG1 = KokkosDeviceView<DGVector<1>>;
    using HostViewDG1 = KokkosHostView<DGVector<1>>;

    // strain and stress components
    using DeviceViewStress = KokkosDeviceView<DGVector<DGstressComp>>;
    using HostViewStress = KokkosHostView<DGVector<DGstressComp>>;
    using ConstDeviceViewStress = ConstKokkosDeviceView<DGVector<DGstressComp>>;

    using DeviceViewAdvect = KokkosDeviceView<DGVector<DGadvection>>;
    using HostViewAdvect = KokkosHostView<DGVector<DGadvection>>;
    using ConstDeviceViewAdvect = ConstKokkosDeviceView<DGVector<DGadvection>>;

    // constant matrices
    static constexpr int NGP = NGP_DG<DGstressComp>;
    using PSIAdvectType = decltype(PSI<DGadvection, NGP>);
    using PSIStressType = decltype(PSI<DGstressComp, NGP>);
    // in gcc13 the signature of updateStressHighOrder is somehow incompatible between declaration
    // and implementation if we use ConstKokkosDeviceView directly
    using PSIAdvectView = ConstKokkosDeviceView<PSIAdvectType>;
    using PSIStressView = ConstKokkosDeviceView<PSIStressType>;

    using EdgeVec = Eigen::Matrix<FloatType, 1, NGP * NGP>;

    // precomputed maps
    using DivMapDevice
        = KokkosDeviceMapView<typename ParametricMomentumMap<CGdegree, DGadvection>::DivMatrix>;
    using GradMapDevice
        = KokkosDeviceMapView<typename ParametricMomentumMap<CGdegree, DGadvection>::GradMatrix>;
    using GaussMapDevice = KokkosDeviceMapView<
        typename ParametricMomentumMap<CGdegree, DGadvection>::GaussMapMatrix>;
    using GaussMapAdvectDevice = KokkosDeviceMapView<
        typename ParametricMomentumMap<CGdegree, DGadvection>::GaussMapAdvectMatrix>;
    using DSSHDevice
        = KokkosDeviceMapView<typename ParametricMomentumMap<CGdegree, DGadvection>::DSSHMatrix>;

    KokkosCGDynamicsKernel(const DynamicsParameters& params);
    // still defaulted but explicitly defined in the source file to allow for pimpl with unique_ptr
    ~KokkosCGDynamicsKernel() override;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    ModelArray getDG0Data(const std::string& name) const override;

    // The host variant is needed in IDynamics:setData where data does not come out of the
    // ModelArrayStore.
    void setData(const std::string& name, const ModelArray& data) override;
    // Use this to directly set the data of kernel's internal device buffers.
    virtual void setData(const std::string& name, const ConstDeviceViewMA& data);
    void setDGArray(const std::string& name, ModelArray::DataType& dgData) override;
    virtual void setDGArray(
        const std::string& name, const KokkosDeviceView<ModelArray::DataType>& dgData);

    void prepareAdvection() override;
    void advectDynamicsFields(FloatType timestep) override;
    DGVector<DGadvection>& advectDGVField(FloatType timestep, DGVector<DGadvection>& field,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity()) override;
    void advectDGVFieldDevice(FloatType timestep, const DeviceViewAdvect& field,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity());

    void updateGradientOfSeaSurfaceHeight();

    // cuda requires these functions to be public but they should only be needed by the concrete
    // dynamics kernel (like protected)
    static void prepareIterationDevice(const DeviceViewCG& cgHDevice, const DeviceViewCG& cgADevice,
        const ConstDeviceViewAdvect& hiceDevice, const ConstDeviceViewAdvect& ciceDevice,
        const Interpolations::KokkosDG2CGInterpolator<CGdegree, DGadvection>& dG2CGInterpolator);

    static void dirichletZero(const DeviceViewCG& v, const ConstDeviceBitset& cgLandMaskDevice);

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
        const ConstDeviceBitset& cgLandMaskDevice, DeviceIndex nx, DeviceIndex ny,
        COORDINATES coordinates);

    static void applyBoundariesDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const ConstDeviceBitset& cgLandMaskDevice, DeviceIndex nx, DeviceIndex ny);

    static void updateIceOceanStressDevice(const DeviceViewCG& uIceOceanStressDevice,
        const DeviceViewCG& vIceOceanStressDevice, const ConstDeviceViewCG& uIceDevice,
        const ConstDeviceViewCG& vIceDevice, const ConstDeviceViewCG& uOceanDevice,
        const ConstDeviceViewCG& vOceanDevice, const DynamicsParameters& params,
        FloatType cosOceanAngle, FloatType sinOceanAngle);

    // functionality from Tools.hpp
    static void computeShearDevice(const DeviceViewAdvect& destDevice,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const PSIStressView& PSIStressDevice,
        const ConstDeviceBitset& landMaskDevice, const GaussMapAdvectDevice& iMJwPSIAdvectDevice);

    static void computeTensorInvariantIDevice(const DeviceViewStress& destDevice,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device);

    static void computeTensorInvariantIIDevice(const DeviceViewStress& destDevice,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device);

protected:
    // currently not used
    void updateMomentum(const TimestepTime& tst) override { }

    // copy data from a ModelArray view to a cG-field on device
    template <typename... Args>
    void mA2CG(const DeviceViewCG& dest,
        const ConstKokkosEigenView<ModelArray::DataType, Args...>& src) const
    {
        Kokkos::deep_copy(tempDataAdvectDevice, 0.0);
        kokkosMA2DG<DGadvection>(tempDataAdvectDevice, src);
        (*dG2CGAdvectInterpolator)(dest, tempDataAdvectDevice);
    }

    // copy data from a dg field on device to a host ModelArray
    ModelArray& dG2MA(ModelArray& ma, const ConstDeviceViewAdvect& dg) const;

    // named fields for setData
    std::unordered_map<std::string, DeviceViewCG> namedCGFields;

    // cG (velocity) components
    DeviceViewCG uDevice;
    HostViewCG uHost;
    DeviceViewCG vDevice;
    HostViewCG vHost;

    DeviceViewCG cgHDevice;
    HostViewCG cgHHost;
    DeviceViewCG cgADevice;
    HostViewCG cgAHost;

    DeviceViewCG xGradSeaSurfaceHeightDevice;
    HostViewCG xGradSeaSurfaceHeightHost;
    DeviceViewCG yGradSeaSurfaceHeightDevice;
    HostViewCG yGradSeaSurfaceHeightHost;
    DeviceViewDG1 seaSurfaceHeightDevice;
    HostViewDG1 seaSurfaceHeightHost;

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

    DeviceViewCG uIceOceanStressDevice;
    HostViewCG uIceOceanStressHost;
    DeviceViewCG vIceOceanStressDevice;
    HostViewCG vIceOceanStressHost;

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

    // vector used to temporary store computed properties like shear for getDG0Data
    DGVector<DGstressComp> tempData;
    DeviceViewStress tempDataDevice;
    HostViewStress tempDataHost;
    DGVector<DGadvection> tempDataAdvect;
    DeviceViewAdvect tempDataAdvectDevice;
    HostViewAdvect tempDataAdvectHost;
    ModelArray tempDataMA;
    DeviceViewMA tempDataMADevice;

    // precomputed parametric map
    DivMapDevice divS1Device;
    DivMapDevice divS2Device;
    DivMapDevice divMDevice;

    GradMapDevice iMgradXDevice;
    GradMapDevice iMgradYDevice;
    GradMapDevice iMMDevice;

    // data that is needed by the child classes implementing stress and momentum
    DeviceViewAdvect hiceDevice;
    // HostViewAdvect hiceHost;
    DeviceViewAdvect ciceDevice;
    // HostViewAdvect ciceHost;
    DeviceViewAdvect hsnowDevice;
    // HostViewAdvect hsnowHost;

    // constant matrices also need to be available on the GPU
    PSIAdvectView PSIAdvectDevice;
    PSIStressView PSIStressDevice;

    // parametric map precomputed transforms
    // todo: refactor into KokkosParametricMap with switch for precomputed / on-the-fly
    ConstDeviceViewCG lumpedCGMassDevice;
    GaussMapDevice iMJwPSIDevice;
    GaussMapAdvectDevice iMJwPSIAdvectDevice;
    ConstDeviceBitset cgLandMaskDevice;

    // held as a pointer because these objects are initialized by their constructors
    std::unique_ptr<KokkosMesh> meshData;
    std::unique_ptr<Interpolations::KokkosCG2DGInterpolator<DGadvection, CGdegree>>
        cG2DGAdvectInterpolator;
    std::unique_ptr<Interpolations::KokkosDG2CGInterpolator<CGdegree, DGadvection>>
        dG2CGAdvectInterpolator;
    std::unique_ptr<KokkosDGTransport<DGadvection>> dGTransportDevice;
    std::unique_ptr<KokkosSlopeLimiter<DGadvection>> slopeLimiterDevice;

    // sea surface height is always computed in first order
    // parts are only initialized if cG2DGAdvectInterpolator has a different order from <1,1>
    std::unique_ptr<Interpolations::KokkosDG2CGInterpolator<1, 1>> dG2CGFirstOrderInterpolator;
    DSSHDevice _dXSSHDevice;
    DSSHDevice _dYSSHDevice;
    ConstDeviceViewCG1 _lumpedCG1MassDevice;
    DeviceViewCG1 _uGradDevice;
    DeviceViewCG1 _vGradDevice;
};

}

#endif // KOKKOSCGDYNAMICSKERNEL_HPP
