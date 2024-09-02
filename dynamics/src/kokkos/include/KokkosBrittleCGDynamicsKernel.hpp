/*!
 * @file KokkosBrittleCGDynamicsKernel.hpp
 * @date August 28, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSBRITTLECGDYNAMICSKERNEL_HPP
#define KOKKOSBRITTLECGDYNAMICSKERNEL_HPP

#include "../include/MEBParameters.hpp"
#include "KokkosCGDynamicsKernel.hpp"

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

    using PSIAdvectView = typename Base::PSIAdvectView;
    using PSIStressView = typename Base::PSIStressView;

    KokkosBrittleCGDynamicsKernel(const MEBParameters& paramsIn)
        : params(paramsIn)
    {
    }

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    // The brittle rheologies use avgU and avgV to do the advection, not u and v, like mEVP
    void prepareAdvection() override { this->dgtransport->prepareAdvection(avgU, avgV); }

    void update(const TimestepTime& tst) override;

    // expose additional fields
    void setData(const std::string& name, const ModelArray& data) override
    {
        if (name == damageName) {
            DGModelArray::ma2dg(data, damage);
        } else {
            CGDynamicsKernel<DGadvection>::setData(name, data);
        }
    }

    ModelArray getDG0Data(const std::string& name) const override
    {

        if (name == damageName) {
            ModelArray data(ModelArray::Type::H);
            return DGModelArray::dg2ma(damage, data);
        } else {
            return CGDynamicsKernel<DGadvection>::getDG0Data(name);
        }
    }

    ModelArray getDGData(const std::string& name) const override
    {
        if (name == damageName) {
            ModelArray data(ModelArray::Type::DG);
            return DGModelArray::dg2ma(damage, data);
        } else {
            return CGDynamicsKernel<DGadvection>::getDGData(name);
        }
    }

    // todo: move this into KokkosBBMDynamics
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
    // for testing only
    using DynamicsKernel<DGadvection, DGstressComp>::nSteps;
    using DynamicsKernel<DGadvection, DGstressComp>::s11;
    using DynamicsKernel<DGadvection, DGstressComp>::s12;
    using DynamicsKernel<DGadvection, DGstressComp>::s22;
    using DynamicsKernel<DGadvection, DGstressComp>::e11;
    using DynamicsKernel<DGadvection, DGstressComp>::e12;
    using DynamicsKernel<DGadvection, DGstressComp>::e22;
    using DynamicsKernel<DGadvection, DGstressComp>::hice;
    using DynamicsKernel<DGadvection, DGstressComp>::cice;
    using DynamicsKernel<DGadvection, DGstressComp>::smesh;
    using DynamicsKernel<DGadvection, DGstressComp>::deltaT;
    using DynamicsKernel<DGadvection, DGstressComp>::stressDivergence;
    using DynamicsKernel<DGadvection, DGstressComp>::applyBoundaries;
    using DynamicsKernel<DGadvection, DGstressComp>::advectionAndLimits;
    using DynamicsKernel<DGadvection, DGstressComp>::dgtransport;

    using CGDynamicsKernel<DGadvection>::u;
    using CGDynamicsKernel<DGadvection>::v;
    using CGDynamicsKernel<DGadvection>::uAtmos;
    using CGDynamicsKernel<DGadvection>::vAtmos;
    using CGDynamicsKernel<DGadvection>::uOcean;
    using CGDynamicsKernel<DGadvection>::vOcean;
    using CGDynamicsKernel<DGadvection>::prepareIteration;
    using CGDynamicsKernel<DGadvection>::projectVelocityToStrain;
    using CGDynamicsKernel<DGadvection>::dirichletZero;
    using CGDynamicsKernel<DGadvection>::cgH;
    using CGDynamicsKernel<DGadvection>::cgA;
    using CGDynamicsKernel<DGadvection>::dStressX;
    using CGDynamicsKernel<DGadvection>::dStressY;
    using CGDynamicsKernel<DGadvection>::pmap;
    void updateMomentum(const TimestepTime& tst) override
    {
#pragma omp parallel for
        for (size_t i = 0; i < u.rows(); ++i) {
            // FIXME dte_over_mass should include snow in the total mass
            const double dteOverMass = deltaT / (params.rho_ice * cgH(i));
            // Memoized initial velocity values
            const double uIce = u(i);
            const double vIce = v(i);

            const double cPrime
                = cgA(i) * params.F_ocean * std::hypot(uOcean(i) - uIce, vOcean(i) - vIce);

            // FIXME grounding term tauB = cBu[i] / std::hypot(uIce, vIce) + u0
            const double tauB = 0.;
            const double alpha = 1 + dteOverMass * (cPrime * cosOceanAngle + tauB);
            /* FIXME latitude needed for spherical cases
             * const double beta = deltaT * params.fc +
             * dteOverMass * cPrime * std::copysign(sinOceanAngle, lat[i]);
             */
            const double beta = deltaT * params.fc + dteOverMass * cPrime * sinOceanAngle;
            const double rDenom = 1 / (SQR(alpha) + SQR(beta));

            // Atmospheric drag
            const double dragAtm = cgA(i) * params.F_atm * std::hypot(uAtmos(i), vAtmos(i));
            const double tauX = dragAtm * uAtmos(i)
                + cPrime * (uOcean(i) * cosOceanAngle - vOcean(i) * sinOceanAngle);
            const double tauY = dragAtm * vAtmos(i)
                + cPrime * (vOcean(i) * cosOceanAngle + uOcean(i) * sinOceanAngle);

            // Stress gradient
            const double gradX = dStressX(i) / pmap->lumpedcgmass(i);
            const double gradY = dStressY(i) / pmap->lumpedcgmass(i);

            u(i) = alpha * uIce + beta * vIce
                + dteOverMass * (alpha * (gradX + tauX) + beta * (gradY + tauY));
            u(i) *= rDenom;

            v(i) = alpha * vIce - beta * uIce
                + dteOverMass * (alpha * (gradY + tauY) + beta * (gradX + tauX));
            v(i) *= rDenom;

            // Calculate the contribution to the average velocity
            avgU(i) += u(i) / nSteps;
            avgV(i) += v(i) / nSteps;
        }
    }

    // these values are only needed on the host because of the advection step
    CGVector<CGdegree> avgU;
    CGVector<CGdegree> avgV;

    DGVector<DGadvection> damage;

    const MEBParameters& params;

    std::unique_ptr<Nextsim::DGTransport<DGstressComp>> stressTransport;

    DeviceViewCG avgUDevice;
    HostViewCG avgUHost;
    DeviceViewCG avgVDevice;
    HostViewCG avgVHost;

    DeviceViewAdvect damageDevice;
    HostViewAdvect damageHost;

    FloatType cosOceanAngle, sinOceanAngle;

    // BBM stress update related
    KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapAdvectMatrix> iMJwPSIAdvectDevice;
    // ParametricMesh::h()
    KokkosDeviceMapView<FloatType> cellSizeDevice;
};

} /* namespace Nextsim */

#endif /* KOKKOSMEVPDYNAMICSKERNEL_HPP */