/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSBRITTLECGDYNAMICSKERNEL_HPP
#define KOKKOSBRITTLECGDYNAMICSKERNEL_HPP

#include "../include/BBMParameters.hpp"
#include "KokkosCGDynamicsKernel.hpp"
#include "KokkosDGTransport.hpp"

namespace Nextsim {

template <int DGadvection>
class KokkosBrittleCGDynamicsKernel : public KokkosCGDynamicsKernel<DGadvection> {
public:
    using Base = KokkosCGDynamicsKernel<DGadvection>;

    using DeviceViewCG = typename Base::DeviceViewCG;
    using HostViewCG = typename Base::HostViewCG;
    using ConstDeviceViewCG = typename Base::ConstDeviceViewCG;

    using DeviceViewStress = typename Base::DeviceViewStress;
    using ConstDeviceViewStress = typename Base::ConstDeviceViewStress;

    using DeviceViewAdvect = typename Base::DeviceViewAdvect;
    using HostViewAdvect = typename Base::HostViewAdvect;
    using ConstDeviceViewAdvect = typename Base::ConstDeviceViewAdvect;

    KokkosBrittleCGDynamicsKernel(const BBMParameters& paramsIn);

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    void advectDynamicsFields(double timestep) override;

    void update(const TimestepTime& tst) override;

    // expose additional fields
    void setData(const std::string& name, const ModelArray& data);
    void setData(const std::string& name, const ConstDeviceViewMA& data);
    // void setDGArray(const std::string& name, ModelArray::DataType& dgData) override;
    using Base::setDGArray;
    void setDGArray(const std::string& name, const DeviceViewMA& dgData) override;

    static void updateMomentumDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const DeviceViewCG& avgUDevice, const DeviceViewCG& avgVDevice,
        const ConstDeviceViewCG& cgHDevice, const ConstDeviceViewCG& cgADevice,
        const ConstDeviceViewCG& uAtmosDevice, const ConstDeviceViewCG& vAtmosDevice,
        const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
        const ConstDeviceViewCG& dStressXDevice, const ConstDeviceViewCG& dStressYDevice,
        const ConstDeviceViewCG& xGradSeaSurfaceHeightDevice,
        const ConstDeviceViewCG& yGradSeaSurfaceHeightDevice,
        const ConstDeviceViewCG& lumpedCGMassDevice, const ConstDeviceBitset& cgLandMaskDevice,
        const FloatType deltaT, const BBMParameters& params, FloatType cosOceanAngle,
        FloatType sinOceanAngle, DeviceIndex nSteps);

protected:
    virtual void updateStressHighOrderDevice(const DeviceViewStress& s11Device,
        const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const ConstDeviceViewAdvect& hiceDevice,
        const ConstDeviceViewAdvect& ciceDevice, const DeviceViewAdvect& damageDevice,
        const FloatType deltaT)
        = 0;

    // Average velocities are only used internally so they don't have to be mirrored on the host.
    // They might still be useful for debugging and outputs.
    CGVector<CGdegree> avgU;
    CGVector<CGdegree> avgV;

    DGVectorHolder<DGadvection> damage;

    const BBMParameters& params;

    std::unique_ptr<KokkosDGTransport<DGstressComp>> stressTransportDevice;
    std::unique_ptr<Interpolations::KokkosCG2DGInterpolator<DGstressComp, CGdegree>>
        cG2DGStressInterpolator;

    DeviceViewCG avgUDevice;
    HostViewCG avgUHost;
    DeviceViewCG avgVDevice;
    HostViewCG avgVHost;

    DeviceViewAdvect damageDevice;
    // HostViewAdvect damageHost;
};

} /* namespace Nextsim */

#endif /* KOKKOSMEVPDYNAMICSKERNEL_HPP */