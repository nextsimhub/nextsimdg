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

    using EdgeVec = Eigen::Matrix<double, 1, NGP * NGP>;

public:
    struct KokkosBuffers {
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
    };

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

        std::tie(buffers.s11Host, buffers.s11Device) = makeKokkosDualView("s11", this->s11);
        std::tie(buffers.s12Host, buffers.s12Device) = makeKokkosDualView("s12", this->s12);
        std::tie(buffers.s22Host, buffers.s22Device) = makeKokkosDualView("s11", this->s22);
        std::tie(buffers.e11Host, buffers.e11Device) = makeKokkosDualView("e11", this->e11);
        std::tie(buffers.e12Host, buffers.e12Device) = makeKokkosDualView("e12", this->e12);
        std::tie(buffers.e22Host, buffers.e22Device) = makeKokkosDualView("e11", this->e22);
        std::tie(buffers.hiceHost, buffers.hiceDevice) = makeKokkosDualView("hice", this->hice);
        std::tie(buffers.ciceHost, buffers.ciceDevice) = makeKokkosDualView("cice", this->cice);

        // does not depend on the data but it can not be allocated in the constructor
        buffers.PSIAdvectDevice
            = makeKokkosDeviceView("PSI<DGadvection, NGP>", PSI<DGadvection, NGP>);
        buffers.PSIStressDevice
            = makeKokkosDeviceView("PSI<DGstress, NGP>", PSI<DGstressDegree, NGP>);

        buffers.iMJwPSIDevice = makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, true);

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

            Kokkos::deep_copy(buffers.s11Device, buffers.s11Host);
            Kokkos::deep_copy(buffers.s12Device, buffers.s12Host);
            Kokkos::deep_copy(buffers.s22Device, buffers.s22Host);
            Kokkos::deep_copy(buffers.e11Device, buffers.e11Host);
            Kokkos::deep_copy(buffers.e12Device, buffers.e12Host);
            Kokkos::deep_copy(buffers.e22Device, buffers.e22Host);
            Kokkos::deep_copy(buffers.hiceDevice, buffers.hiceHost);
            Kokkos::deep_copy(buffers.ciceDevice, buffers.ciceHost);
            std::array<DGVector<DGstressDegree>, N_TENSOR_ELEMENTS> stressTemp { this->s11,
                this->s12, this->s22 };
            std::array<std::reference_wrapper<DGVector<DGstressDegree>>, N_TENSOR_ELEMENTS> stress
                = { stressTemp[0], stressTemp[1],
                      stressTemp[2] };
            stressStep.stressUpdateHighOrder(params, *this->smesh, stress,
                { this->e11, this->e12, this->e22 }, this->hice, this->cice, this->deltaT);

            stressUpdateHighOrder(buffers, params, alpha);

            Kokkos::deep_copy(buffers.s11Host, buffers.s11Device);
            Kokkos::deep_copy(buffers.s12Host, buffers.s12Device);
            Kokkos::deep_copy(buffers.s22Host, buffers.s22Device);

            std::cout << "ref: " << stressTemp[0].norm() << ", gpu: " << this->s11.norm()
                      << ", dif: " << (stressTemp[0] - this->s11).norm() << "\n";

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

    static void stressUpdateHighOrder(
        const KokkosBuffers& _buffers, const VPParameters& _params, double _alpha)
    {
        const int n = _buffers.s11Device.extent(0);
        Kokkos::parallel_for(
            "stressUpdateHighOrder", n, KOKKOS_LAMBDA(const int i) {
                auto s11 = makeEigenMap(_buffers.s11Device);
                auto s12 = makeEigenMap(_buffers.s12Device);
                auto s22 = makeEigenMap(_buffers.s22Device);
                auto e11 = makeEigenMap(_buffers.e11Device);
                auto e12 = makeEigenMap(_buffers.e12Device);
                auto e22 = makeEigenMap(_buffers.e22Device);
                auto hice = makeEigenMap(_buffers.hiceDevice);
                auto cice = makeEigenMap(_buffers.ciceDevice);

                auto PSIAdvect = makeEigenMap(_buffers.PSIAdvectDevice);
                auto PSIStress = makeEigenMap(_buffers.PSIStressDevice);
                /*           auto S11 = AlignedMap<DGVec<DGstress>>(
                              S11Device.data(), nElem, DGVec<DGstress>::ColsAtCompileTime);
                          auto S12 = AlignedMap<DGVec<DGstress>>(
                              S12Device.data(), nElem, DGVec<DGstress>::ColsAtCompileTime);
                          auto S22 = AlignedMap<DGVec<DGstress>>(
                              S22Device.data(), nElem, DGVec<DGstress>::ColsAtCompileTime);

                          auto E11 = AlignedMap<const DGVec<DGstress>>(
                              E11Device.data(), nElem, DGVec<DGstress>::ColsAtCompileTime);
                          auto E12 = AlignedMap<const DGVec<DGstress>>(
                              E12Device.data(), nElem, DGVec<DGstress>::ColsAtCompileTime);
                          auto E22 = AlignedMap<const DGVec<DGstress>>(
                              E22Device.data(), nElem, DGVec<DGstress>::ColsAtCompileTime);

                          auto H = AlignedMap<const DGVec<DGadvection>>(
                              Hdevice.data(), nElem, DGVec<DGadvection>::ColsAtCompileTime);
                          auto A = AlignedMap<const DGVec<DGadvection>>(
                              Adevice.data(), nElem, DGVec<DGadvection>::ColsAtCompileTime);*/

                auto h_gauss = (hice.row(i) * PSIAdvect).array().max(0.0).matrix();
                auto a_gauss = (cice.row(i) * PSIAdvect).array().max(0.0).min(1.0).matrix();

                EdgeVec P
                    = (_params.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp())
                          .matrix();
                const EdgeVec e11_gauss = e11.row(i) * PSIStress;
                const EdgeVec e12_gauss = e12.row(i) * PSIStress;
                const EdgeVec e22_gauss = e22.row(i) * PSIStress;

                const auto DELTA = (_params.DeltaMin * _params.DeltaMin
                    + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                    + 1.50 * e11_gauss.array() * e22_gauss.array() + e12_gauss.array().square())
                                       .sqrt()
                                       .matrix();

                const auto& map = _buffers.iMJwPSIDevice[i];

                const FloatType alphaInv = 1.0 / _alpha;
                const FloatType fac = 1.0 - alphaInv;
                const EdgeVec PDelta = P.array() / DELTA.array();
                s11.row(i) = fac * s11.row(i)
                    + (map
                        * (alphaInv
                            * (PDelta.array()
                                    * ((5.0 / 8.0) * e11_gauss.array()
                                        + (3.0 / 8.0) * e22_gauss.array())
                                - 0.5 * P.array())
                                  .matrix()
                                  .transpose()))
                          .transpose();
                s22.row(i) = fac * s22.row(i)
                    + (map
                        * (alphaInv
                            * (PDelta.array()
                                    * ((5.0 / 8.0) * e22_gauss.array()
                                        + (3.0 / 8.0) * e11_gauss.array())
                                - 0.5 * P.array())
                                  .matrix()
                                  .transpose()))
                          .transpose();
                s12.row(i) = fac * s12.row(i)
                    + (map
                        * (alphaInv
                            * (PDelta.array() * (1.0 / 4.0) * e12_gauss.array())
                                  .matrix()
                                  .transpose()))
                          .transpose();
            });
    }

private:
    MEVPStressUpdateStep<DGadvection, DGstressDegree, CGdegree> stressStep;
    KokkosBuffers buffers;
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
