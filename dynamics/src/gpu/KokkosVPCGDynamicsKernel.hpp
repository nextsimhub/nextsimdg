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
    //        , PSIAdvectDevice(makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection,
    //        NGP>)) , PSIStressDevice(makeKokkosDeviceView("PSI<DGstress, NGP>",
    //        PSI<DGstressDegree, NGP>))
    {
    }

    KokkosVPCGDynamicsKernel(const KokkosVPCGDynamicsKernel<DGadvection>&) = delete;
    KokkosVPCGDynamicsKernel(KokkosVPCGDynamicsKernel<DGadvection>&&) = delete;

    KokkosVPCGDynamicsKernel& operator=(const KokkosVPCGDynamicsKernel<DGadvection>&) = delete;
    KokkosVPCGDynamicsKernel& operator=(KokkosVPCGDynamicsKernel<DGadvection>&&) = delete;

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override
    {
        CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

        std::tie(s11Host, s11Device) = makeKokkosDualView("s11", this->s11);
        std::tie(s12Host, s12Device) = makeKokkosDualView("s12", this->s12);
        std::tie(s22Host, s22Device) = makeKokkosDualView("s11", this->s22);
        std::tie(e11Host, e11Device) = makeKokkosDualView("e11", this->e11);
        std::tie(e12Host, e12Device) = makeKokkosDualView("e12", this->e12);
        std::tie(e22Host, e22Device) = makeKokkosDualView("e11", this->e22);
        std::tie(hiceHost, hiceDevice) = makeKokkosDualView("hice", this->hice);
        std::tie(ciceHost, ciceDevice) = makeKokkosDualView("cice", this->cice);

        // does not depend on the data but it can not be allocated in the constructor
        PSIAdvectDevice = makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection, NGP>);
        PSIStressDevice = makeKokkosDeviceView("PSI<DGstress, NGP>", PSI<DGstressDegree, NGP>);

        iMJwPSIDevice = makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, true);

        stressStep.setPMap(this->pmap);
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

            Kokkos::deep_copy(s11Device, s11Host);
            Kokkos::deep_copy(s12Device, s12Host);
            Kokkos::deep_copy(s22Device, s22Host);
            Kokkos::deep_copy(e11Device, e11Host);
            Kokkos::deep_copy(e12Device, e12Host);
            Kokkos::deep_copy(e22Device, e22Host);
            Kokkos::deep_copy(hiceDevice, hiceHost);
            Kokkos::deep_copy(ciceDevice, ciceHost);
            std::array<DGVector<DGstressDegree>, N_TENSOR_ELEMENTS> stressTemp { this->s11,
                this->s12, this->s22 };
            std::array<std::reference_wrapper<DGVector<DGstressDegree>>, N_TENSOR_ELEMENTS> stress
                = { stressTemp[0], stressTemp[1],
                      stressTemp[2] }; // Call the step function on the StressUpdateStep class
            // Call the step function on the StressUpdateStep class
            stressStep.stressUpdateHighOrder(params, *this->smesh, stress,
                { this->e11, this->e12, this->e22 }, this->hice, this->cice, this->deltaT);

            Kokkos::deep_copy(s11Host, s11Device);
            Kokkos::deep_copy(s12Host, s12Device);
            Kokkos::deep_copy(s22Host, s22Device);
            Kokkos::deep_copy(e11Host, e11Device);
            Kokkos::deep_copy(e12Host, e12Device);
            Kokkos::deep_copy(e22Host, e22Device);
            Kokkos::deep_copy(hiceHost, hiceDevice);
            Kokkos::deep_copy(ciceHost, ciceDevice);

            double stressScale
                = 1.0; // 2nd-order Stress term has different scaling with the EarthRadius
            if (this->smesh->CoordinateSystem == Nextsim::SPHERICAL)
                stressScale = 1.0 / Nextsim::EarthRadius / Nextsim::EarthRadius;

            this->stressDivergence(); // Compute divergence of stress tensor

            updateMomentum(tst);

            this->applyBoundaries();
        }
        // Finally, do the base class update
        DynamicsKernel<DGadvection, DGstressDegree>::update(tst);
    }

private:
    MEVPStressUpdateStep<DGadvection, DGstressDegree, CGdegree> stressStep;
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
