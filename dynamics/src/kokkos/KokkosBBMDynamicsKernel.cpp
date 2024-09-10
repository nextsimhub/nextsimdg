/*!
 * @file KokkosBrittleCGDynamicsKernel.cpp
 * @date August 28, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosBBMDynamicsKernel.hpp"

namespace Nextsim {

template <int DGadvection>
KokkosBBMDynamicsKernel<DGadvection>::KokkosBBMDynamicsKernel(const MEBParameters& paramsIn)
    : KokkosBrittleCGDynamicsKernel<DGadvection>(paramsIn)
{
}

template <int DGadvection>
void KokkosBBMDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    KokkosBrittleCGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    iMJwPSIAdvectDevice = makeKokkosDeviceViewMap("iMJwPSIAdvect", this->pmap->iMJwPSIAdvect, true);
    std::vector<FloatType> cellSize(this->smesh->nelements);
    for (size_t i = 0; i < this->smesh->nelements; ++i) {
        cellSize[i] = this->smesh->h(i);
    }
    cellSizeDevice = makeKokkosDeviceViewMap("cellSize", cellSize, true);
}

template <int DGadvection>
void KokkosBBMDynamicsKernel<DGadvection>::updateStressHighOrderDevice(
    const DeviceViewStress& s11Device, const DeviceViewStress& s12Device,
    const DeviceViewStress& s22Device, const ConstDeviceViewStress& e11Device,
    const ConstDeviceViewStress& e12Device, const ConstDeviceViewStress& e22Device,
    const ConstDeviceViewAdvect& hiceDevice, const ConstDeviceViewAdvect& ciceDevice,
    const DeviceViewAdvect& damageDevice, const FloatType deltaT)
{
    updateStressHighOrderDevice(s11Device, s12Device, s22Device, e11Device, e12Device, e22Device,
        this->PSIAdvectDevice, this->PSIStressDevice, hiceDevice, ciceDevice, damageDevice,
        this->iMJwPSIDevice, this->iMJwPSIAdvectDevice, this->cellSizeDevice, deltaT, this->params);
}

template <int DGadvection>
void KokkosBBMDynamicsKernel<DGadvection>::updateStressHighOrderDevice(
    const DeviceViewStress& s11Device, const DeviceViewStress& s12Device,
    const DeviceViewStress& s22Device, const ConstDeviceViewStress& e11Device,
    const ConstDeviceViewStress& e12Device, const ConstDeviceViewStress& e22Device,
    const PSIAdvectView& PSIAdvectDevice, const PSIStressView& PSIStressDevice,
    const ConstDeviceViewAdvect& hiceDevice, const ConstDeviceViewAdvect& ciceDevice,
    const DeviceViewAdvect& damageDevice,
    const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapMatrix>& iMJwPSIDevice,
    const KokkosDeviceMapView<ParametricMomentumMap<CGdegree>::GaussMapAdvectMatrix>&
        iMJwPSIAdvectDevice,
    const KokkosDeviceMapView<FloatType>& cellSizeDevice, const FloatType deltaT,
    const MEBParameters& params)
{
    constexpr int NGP = KokkosCGDynamicsKernel<DGadvection>::NGP;
    using EdgeVec = Eigen::Matrix<FloatType, 1, NGP * NGP>;

    Kokkos::parallel_for(
        "updateStressHighOrder", s11Device.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto s11 = makeEigenMap(s11Device);
            auto s12 = makeEigenMap(s12Device);
            auto s22 = makeEigenMap(s22Device);
            auto e11 = makeEigenMap(e11Device);
            auto e12 = makeEigenMap(e12Device);
            auto e22 = makeEigenMap(e22Device);
            auto hice = makeEigenMap(hiceDevice);
            auto cice = makeEigenMap(ciceDevice);
            auto damage = makeEigenMap(damageDevice);

            const auto PSIAdvect = makeEigenMap(PSIAdvectDevice);
            const auto PSIStress = makeEigenMap(PSIStressDevice);

            const EdgeVec hGauss = (hice.row(i) * PSIAdvect).array().max(0.0).matrix();
            const EdgeVec aGauss = (cice.row(i) * PSIAdvect).array().max(0.0).min(1.0).matrix();
            EdgeVec dGauss = (damage.row(i) * PSIAdvect).array().max(1e-12).min(1.0).matrix();

            const EdgeVec e11Gauss = e11.row(i) * PSIStress;
            const EdgeVec e12Gauss = e12.row(i) * PSIStress;
            const EdgeVec e22Gauss = e22.row(i) * PSIStress;

            EdgeVec s11Gauss = s11.row(i) * PSIStress;
            EdgeVec s12Gauss = s12.row(i) * PSIStress;
            EdgeVec s22Gauss = s22.row(i) * PSIStress;

            //! Current normal stress for the evaluation of tildeP (Eqn. 1)
            EdgeVec sigma_n = 0.5 * (s11Gauss.array() + s22Gauss.array());

            //! exp(-C(1-A))
            const EdgeVec expC = (params.compaction_param * (1.0 - aGauss.array())).exp().array();

            // Eqn. 25
            const auto powalphaexpC
                = (dGauss.array() * expC.array()).pow(params.exponent_relaxation_sigma - 1);
            const EdgeVec time_viscous = params.undamaged_time_relaxation_sigma * powalphaexpC;

            //! BBM  Computing tildeP according to (Eqn. 7b and Eqn. 8)
            // (Eqn. 8)
            const auto Pmax
                = params.P0 * hGauss.array().pow(params.exponent_compression_factor) * expC.array();

            // (Eqn. 7b) Prepare tildeP
            // tildeP must be capped at 1 to get an elastic response
            // (Eqn. 7b) Select case based on sigma_n
            const auto tildeP
                = (sigma_n.array() < 0.0)
                      .select((-Pmax.array() / sigma_n.array()).min(1.0).matrix(), 0.);

            // multiplicator
            const EdgeVec multiplicator
                = time_viscous.array() / (time_viscous.array() + (1. - tildeP.array()) * deltaT);

            //! Eqn. 9
            const EdgeVec elasticity
                = hGauss.array() * params.young * dGauss.array() * expC.array();

            // Eqn. 12: first factor on RHS
            /* Stiffness matrix
             * / (K:e)11 \       1     /  1  nu    0  \ / e11 \
             * | (K:e)22 |  =  ------- | nu   1    0  | | e22 |
             * \ (K:e)12 /    1 - nu^2 \  0   0  1-nu / \ e12 /
             */

            const EdgeVec Dunit_factor
                = deltaT * elasticity.array() / (1. - (params.nu0 * params.nu0));

            s11Gauss.array()
                += Dunit_factor.array() * (e11Gauss.array() + params.nu0 * e22Gauss.array());
            s22Gauss.array()
                += Dunit_factor.array() * (params.nu0 * e11Gauss.array() + e22Gauss.array());
            s12Gauss.array() += Dunit_factor.array() * e12Gauss.array() * (1. - params.nu0);

            // //! Implicit part of RHS (Eqn. 33)
            s11Gauss.array() *= multiplicator.array();
            s22Gauss.array() *= multiplicator.array();
            s12Gauss.array() *= multiplicator.array();

            sigma_n = 0.5 * (s11Gauss.array() + s22Gauss.array());
            const EdgeVec tau = (0.25 * (s11Gauss.array() - s22Gauss.array()).square()
                + s12Gauss.array().square())
                                    .sqrt();

            const FloatType scale_coef = Kokkos::sqrt(0.1 / cellSizeDevice(i));

            //! Eqn. 22
            const auto cohesion = params.C_lab * scale_coef * hGauss.array();
            //! Eqn. 30
            const EdgeVec compr_strength = params.compr_strength * scale_coef * hGauss.array();

            // Mohr-Coulomb failure using Mssrs. Plante & Tremblay's formulation
            // sigma_s + tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
            EdgeVec dcrit
                = (tau.array() + params.tan_phi * sigma_n.array() > 0.)
                      .select(
                          cohesion.array() / (tau.array() + params.tan_phi * sigma_n.array()), 1.);

            // Compressive failure using Mssrs. Plante & Tremblay's formulation
            dcrit = (sigma_n.array() < -compr_strength.array())
                        .select(-compr_strength.array() / sigma_n.array(), dcrit);

            // Only damage when we're outside
            dcrit = dcrit.array().min(1.0);

            // Eqn. 29
            const auto td = cellSizeDevice(i)
                * Kokkos::sqrt(2. * (1. + params.nu0) * params.rho_ice) / elasticity.array().sqrt();

            // Update damage
            const EdgeVec relaxFac = (1. - dcrit.array()) * deltaT / td.array();
            dGauss.array() -= dGauss.array() * relaxFac.array();

            // Relax stress in Gauss points
            s11Gauss.array() -= s11Gauss.array() * relaxFac.array();
            s12Gauss.array() -= s12Gauss.array() * relaxFac.array();
            s22Gauss.array() -= s22Gauss.array() * relaxFac.array();

            // INTEGRATION OF STRESS AND DAMAGE
            // get the inverse of the mass matrix scaled with the test-functions in the gauss
            // points, with the gauss weights and with J. This is a 8 x 9 matrix
            const auto iMJwPSI = iMJwPSIDevice[i];
            s11.row(i) = iMJwPSI * s11Gauss.matrix().transpose();
            s12.row(i) = iMJwPSI * s12Gauss.matrix().transpose();
            s22.row(i) = iMJwPSI * s22Gauss.matrix().transpose();

            damage.row(i) = iMJwPSIAdvectDevice[i] * dGauss.matrix().transpose();
        });
}

// because ParametricMomentumMap<CGdegree>::iMJwPSIAdvect does not properly depend on DGadvection we
// can only build this version since the switch is implemented in compile-time we dont really need
// the other versions anyway
template class KokkosBBMDynamicsKernel<DGCOMP>;

}