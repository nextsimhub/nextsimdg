/*!
 * @file KokkosBrittleCGDynamicsKernel.hpp
 * @date August 28, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSBRITTLECGDYNAMICSKERNEL_HPP
#define KOKKOSBRITTLECGDYNAMICSKERNEL_HPP

#include "../include/MEBParameters.hpp"
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

    KokkosBrittleCGDynamicsKernel(const MEBParameters& paramsIn);

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    // The brittle rheologies use avgU and avgV to do the advection, not u and v, like mEVP
    void prepareAdvection() override { this->dgtransport->prepareAdvection(avgU, avgV); }

    void update(const TimestepTime& tst) override;

    // expose additional fields
    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) const override;
    ModelArray getDGData(const std::string& name) const override;

    static void updateMomentumDevice(const DeviceViewCG& uDevice, const DeviceViewCG& vDevice,
        const DeviceViewCG& avgUDevice, const DeviceViewCG& avgVDevice,
        const ConstDeviceViewCG& cgHDevice, const ConstDeviceViewCG& cgADevice,
        const ConstDeviceViewCG& uAtmosDevice, const ConstDeviceViewCG& vAtmosDevice,
        const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
        const ConstDeviceViewCG& dStressXDevice, const ConstDeviceViewCG& dStressYDevice,
        const ConstDeviceViewCG& lumpedCGMassDevice, const FloatType deltaT,
        const MEBParameters& params, FloatType cosOceanAngle, FloatType sinOceanAngle,
        DeviceIndex nSteps);

protected:
    virtual void updateStressHighOrderDevice(const DeviceViewStress& s11Device,
        const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const ConstDeviceViewAdvect& hiceDevice,
        const ConstDeviceViewAdvect& ciceDevice, const DeviceViewAdvect& damageDevice,
        const FloatType deltaT)
        = 0;

    // these values are only needed on the host because of the advection step
    CGVector<CGdegree> avgU;
    CGVector<CGdegree> avgV;

    DGVector<DGadvection> damage;

    const MEBParameters& params;

    std::unique_ptr<DGTransport<DGstressComp>> stressTransport;
    std::unique_ptr<KokkosDGTransport<DGstressComp>> stressTransportDevice;
    std::unique_ptr<Interpolations::KokkosCG2DGInterpolator<DGstressComp, CGdegree>>
        cG2DGStressInterpolator;

    DeviceViewCG avgUDevice;
    HostViewCG avgUHost;
    DeviceViewCG avgVDevice;
    HostViewCG avgVHost;

    DeviceViewAdvect damageDevice;
    HostViewAdvect damageHost;

    FloatType cosOceanAngle, sinOceanAngle;
};

} /* namespace Nextsim */

#endif /* KOKKOSMEVPDYNAMICSKERNEL_HPP */