/*!
 * @file KokkosBBMDynamicsKernel.hpp
 * @date August 28, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSBBMDYNAMICSKERNEL_HPP
#define KOKKOSBBMDYNAMICSKERNEL_HPP

#include "KokkosBrittleCGDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection>
class KokkosBBMDynamicsKernel : public KokkosBrittleCGDynamicsKernel<DGadvection> {
public:
    using Base = KokkosBrittleCGDynamicsKernel<DGadvection>;

    using DeviceViewStress = typename Base::DeviceViewStress;
    using ConstDeviceViewStress = typename Base::ConstDeviceViewStress;
    using DeviceViewAdvect = typename Base::DeviceViewAdvect;
    using ConstDeviceViewAdvect = typename Base::ConstDeviceViewAdvect;

    using PSIAdvectView = typename KokkosCGDynamicsKernel<DGadvection>::PSIAdvectView;
    using PSIStressView = typename KokkosCGDynamicsKernel<DGadvection>::PSIStressView;

    KokkosBBMDynamicsKernel(const MEBParameters& paramsIn);

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    static void updateStressHighOrderDevice(const DeviceViewStress& s11Device,
        const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const PSIAdvectView& PSIAdvectDevice,
        const PSIStressView& PSIStressDevice, const ConstDeviceViewAdvect& hiceDevice,
        const ConstDeviceViewAdvect& ciceDevice, const DeviceViewAdvect& damageDevice,
        const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix>& iMJwPSIDevice,
        const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapAdvectMatrix>&
            iMJwPSIAdvectDevice,
        const KokkosDeviceMapView<FloatType>& cellSizeDevice, const FloatType deltaT,
        const MEBParameters& params);

protected:
    void updateStressHighOrderDevice(const DeviceViewStress& s11Device,
        const DeviceViewStress& s12Device, const DeviceViewStress& s22Device,
        const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
        const ConstDeviceViewStress& e22Device, const ConstDeviceViewAdvect& hiceDevice,
        const ConstDeviceViewAdvect& ciceDevice, const DeviceViewAdvect& damageDevice,
        const FloatType deltaT) override;

private:
    KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapAdvectMatrix> iMJwPSIAdvectDevice;
    // ParametricMesh::h()
    KokkosDeviceMapView<FloatType> cellSizeDevice;
};

} /* namespace Nextsim */

#endif /* KOKKOSBBMDYNAMICSKERNEL_HPP */