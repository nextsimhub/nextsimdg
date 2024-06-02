/*!
 * @file KokkosVPCGDynamicsKernel.hpp
 *
 * @date Feb 2, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSVPCGDYNAMICSKERNEL_HPP
#define KOKKOSVPCGDYNAMICSKERNEL_HPP

#include "../include/CGDynamicsKernel.hpp"
#include "../include/VPParameters.hpp"
#include "KokkosUtils.hpp"

#include <Kokkos_Bitset.hpp>

// only for testing
#include "../include/MEVPStressUpdateStep.hpp"

namespace Nextsim {

template <int DG> constexpr int NGP_DG = ((DG == 8) || (DG == 6)) ? 3 : (DG == 3 ? 2 : -1);

// The VP pseudo-timestepping momentum equation solver for CG velocities
template <int DGadvection> class KokkosVPCGDynamicsKernel : public CGDynamicsKernel<DGadvection> {
private:
    static constexpr int NGP = NGP_DG<DGstressDegree>;

    using EdgeVec = Eigen::Matrix<double, 1, NGP * NGP>;

public:
    struct KokkosBuffers {
        // velocity components
        using DeviceViewCG = KokkosDeviceView<CGVector<CGdegree>>;
        using HostViewCG = KokkosHostView<CGVector<CGdegree>>;
        DeviceViewCG uDevice;
        HostViewCG uHost;
        DeviceViewCG vDevice;
        HostViewCG vHost;

        DeviceViewCG dStressXDevice;
        HostViewCG dStressXHost;
        DeviceViewCG dStressYDevice;
        HostViewCG dStressYHost;

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

        // parametric map precomputed transforms
        // todo: refactor into KokkosParametricMap with switch for precomputed / on-the-fly
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix> divS1Device;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix> divS2Device;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::DivMatrix> divMDevice;
        
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix> iMgradXDevice;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix> iMgradYDevice;
        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GradMatrix> iMMDevice;

        KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix> iMJwPSIDevice;

        // mesh related
        std::array<KokkosDeviceMapView<size_t>,4> dirichletDevice;
        Kokkos::ConstBitset<Kokkos::DefaultExecutionSpace> landMaskDevice;
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

    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void update(const TimestepTime& tst) override;

    // todo: move kokkos stuff out of class into extra namespace?
    // cuda requires these functions to be public
    static void projVelocityToStrain(
        const KokkosBuffers& _buffers, int nx, int ny, COORDINATES coordinates);
    static void stressUpdateHighOrder(
        const KokkosBuffers& _buffers, const VPParameters& _params, double _alpha);
    static void computeStressDivergence(const KokkosBuffers& _buffers, int nx, int ny, COORDINATES coordinates);

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
