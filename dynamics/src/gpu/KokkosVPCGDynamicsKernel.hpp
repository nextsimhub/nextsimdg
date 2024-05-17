/*!
 * @file KokkosVPCGDynamicsKernel.hpp
 *
 * @date Feb 2, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSVPCGDYNAMICSKERNEL_HPP
#define KOKKOSVPCGDYNAMICSKERNEL_HPP

#include "../include/CGDynamicsKernel.hpp"
#include "KokkosUtils.hpp"

namespace Nextsim {

template <int DG> constexpr int NGP_DG = ((DG == 8) || (DG == 6)) ? 3 : (DG == 3 ? 2 : -1);

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection> class KokkosVPCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
private:
    static constexpr int NGP = NGP_DG<DGstressDegree>;

    // strain and stress components
    using DeviceViewStress = KokkosDeviceView<DGVector<DGstressDegree>>;
    using HostViewStress = KokkosHostView<DGVector<DGstressDegree>>;
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

    // constant matrices also need to be available on the GPU
    using PSIAdvectType = decltype(PSI<DGadvection, NGP>);
    using PSIStressType = decltype(PSI<DGstressDegree, NGP>);
    ConstKokkosDeviceView<PSIAdvectType> PSIAdvectDevice;
    ConstKokkosDeviceView<PSIStressType> PSIStressDevice;

    using IMJwPSIType = Eigen::Matrix<FloatType, DGstressDegree, NGP * NGP>;
    KokkosDeviceMapView<IMJwPSIType> iMJwPSIDevice;

public:
    KokkosVPCGDynamicsKernel(const VPParameters& paramsIn)
        : CGDynamicsKernel<DGadvection>()
        , params(paramsIn)
        , s11Device(makeKokkosDeviceView("s11", this->s11))
        , s11Host(makeKokkosHostView(this->s11))
        , s12Device(makeKokkosDeviceView("s12", this->s12))
        , s12Host(makeKokkosHostView(this->s12))
        , s22Device(makeKokkosDeviceView("s22", this->s22))
        , s22Host(makeKokkosHostView(this->s22))
        , e11Device(makeKokkosDeviceView("e11", this->e11))
        , e11Host(makeKokkosHostView(this->e11))
        , e12Device(makeKokkosDeviceView("e12", this->e12))
        , e12Host(makeKokkosHostView(this->e12))
        , e22Device(makeKokkosDeviceView("e22", this->e22))
        , e22Host(makeKokkosHostView(this->e22))
        , hiceDevice(makeKokkosDeviceView("hice", this->hice))
        , hiceHost(makeKokkosHostView( this->hice))
        , ciceDevice(makeKokkosDeviceView("cice", this->cice))
        , ciceHost(makeKokkosHostView( this->cice))
        , PSIAdvectDevice(makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection, NGP>))
        , PSIStressDevice(makeKokkosDeviceView("PSI<DGstress, NGP>", PSI<DGstressDegree, NGP>))
        , iMJwPSIDevice(makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, true))
    {
    }

    void update(const TimestepTime& tst) override
    {
        // Let DynamicsKernel handle the advection step
        DynamicsKernel<DGadvection, DGstressDegree>::advectionAndLimits(tst);

        this->prepareIteration({ { hiceName, this->hice }, { ciceName, this->cice } });

        u0 = this->u;
        v0 = this->v;

        // The critical timestep for the VP solver is the advection timestep
        this->deltaT = tst.step.seconds();

        for (size_t mevpstep = 0; mevpstep < this->nSteps; ++mevpstep) {
            this->projectVelocityToStrain();

            // Call the step function on the StressUpdateStep class

            double stressScale
                = 1.0; // 2nd-order Stress term has different scaling with the EarthRadius
            if (this->smesh->CoordinateSystem == Nextsim::SPHERICAL)
                stressScale = 1.0 / Nextsim::EarthRadius / Nextsim::EarthRadius;

            this->stressDivergence(stressScale); // Compute divergence of stress tensor

            updateMomentum(tst);

            this->applyBoundaries();
        }
        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressDegree>::update(tst);
    }

private:
    const VPParameters& params;
    double alpha = 1500.;
    double beta = 1500.;

    // Step-initial ice velocity
    CGVector<CGdegree> u0;
    CGVector<CGdegree> v0;

    void updateMomentum(const TimestepTime& tst) override
    {

        // Update the velocity
        double SC = 1.0; ///(1.0-pow(1.0+1.0/beta,-1.0*nSteps));

        //      update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
        for (int i = 0; i < this->u.rows(); ++i) {
            auto uOcnRel = this->uOcean(i)
                - this->u(i); // note the reversed sign compared to the v component
            auto vOcnRel = this->v(i) - this->vOcean(i);
            double absatm = sqrt(SQR(this->uAtmos(i)) + SQR(this->vAtmos(i)));
            double absocn = sqrt(
                SQR(uOcnRel) + SQR(vOcnRel)); // note that the sign of uOcnRel is irrelevant here

            this->u(i) = (1.0
                / (params.rho_ice * this->cgH(i) / this->deltaT * (1.0 + beta) // implicit parts
                    + this->cgA(i) * params.F_ocean * absocn) // implicit parts
                * (params.rho_ice * this->cgH(i) / this->deltaT
                        * (beta * this->u(i) + u0(i)) // pseudo-timestepping
                    + this->cgA(i)
                        * (params.F_atm * absatm * this->uAtmos(i) + // atm forcing
                            params.F_ocean * absocn * SC * this->uOcean(i)) // ocean forcing
                    + params.rho_ice * this->cgH(i) * params.fc * vOcnRel // cor + surface
                    + this->dStressX(i) / this->pmap->lumpedcgmass(i)));
            this->v(i) = (1.0
                / (params.rho_ice * this->cgH(i) / this->deltaT * (1.0 + beta) // implicit parts
                    + this->cgA(i) * params.F_ocean * absocn) // implicit parts
                * (params.rho_ice * this->cgH(i) / this->deltaT
                        * (beta * this->v(i) + v0(i)) // pseudo-timestepping
                    + this->cgA(i)
                        * (params.F_atm * absatm * this->vAtmos(i) + // atm forcing
                            params.F_ocean * absocn * SC * this->vOcean(i)) // ocean forcing
                    + params.rho_ice * this->cgH(i) * params.fc
                        * uOcnRel // here the reversed sign of uOcnRel is used
                    + this->dStressY(i) / this->pmap->lumpedcgmass(i)));
        }
    }
};

} /* namespace Nextsim */

#endif /* KOKKOSVPCGDYNAMICSKERNEL_HPP */
