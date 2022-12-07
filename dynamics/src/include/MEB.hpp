/*!
 * @file MEB.hpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MEB_HPP
#define __MEB_HPP

#include "MEBParameters.hpp"
#include "ParametricTools.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

// Benchmark LKF testcase from [Mehlmann et al., 2022]
namespace RefScale {
constexpr double T = 2.0 * 24 * 60. * 60.; //!< Time horizon 2 days
constexpr double L = 512000.0; //!< Size of domain !!!
constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind
}

// Benchmark compression testcase from [Plante Tremblay, 2021]
namespace RefScaleCompression {
using namespace RefScale;

constexpr double T = 4 * 60. * 60.; //!< Time horizon 2 hours
constexpr double L = 60000.0; //!< Size of domain
constexpr double vmax_ocean = 0.0; //!< Maximum velocity of ocean
constexpr double vmax_atm = 20; //!< Max. vel. of wind

}

namespace Nextsim {

/*!
 * This namespace collects the routines required for the MEB solver
 */
namespace MEB {

    inline constexpr double SQR(double x) { return x * x; }

    /*!
     * @brief Calculate Stresses for the current time step and update damage.
     *
     * @details BBM model, numbering of equations in comments according to:
     *          https://doi.org/10.1029/2021MS002685
     *
     * @tparam CG velocity degree
     * @tparam DGs Stress adn Strain DG degree
     * @tparam DGa Advection H, A, D DG degree
     * @param mesh mesh
     * @param S11 Stress component 11
     * @param S12 Stress component 12
     * @param S22 Stress component 22
     * @param E11 Strain component 11
     * @param E12 Strain component 12
     * @param E22 Strain component 22
     * @param H ice height
     * @param A ice concentation
     * @param D damage
     * @param dt_mom timestep for momentum subcycle
     */
    template <int CG, int DGs, int DGa>
    void StressUpdateHighOrder(const MEBParameters& params,
        const ParametricMesh& smesh, DGVector<DGs>& S11, DGVector<DGs>& S12,
        DGVector<DGs>& S22, const DGVector<DGs>& E11, const DGVector<DGs>& E12,
        const DGVector<DGs>& E22, const DGVector<DGa>& H,
        const DGVector<DGa>& A, DGVector<DGa>& D,
        const double dt_mom)
    {
#define NGP 3
        //! Stress and Damage Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            //! Evaluate values in Gauss points (3 point Gauss rule in 2d => 9 points)
            const Eigen::Matrix<double, 1, NGP* NGP> h_gauss = (H.row(i) * PSI<DGa, NGP>).array().max(0.0).matrix();
            const Eigen::Matrix<double, 1, NGP* NGP> a_gauss = (A.row(i) * PSI<DGa, NGP>).array().max(0.0).min(1.0).matrix();
            Eigen::Matrix<double, 1, NGP* NGP> d_gauss = (D.row(i) * PSI<DGa, NGP>).array().max(0.0).min(1.0-1e-12).matrix();

            const Eigen::Matrix<double, 1, NGP* NGP> e11_gauss = E11.row(i) * PSI<DGs, NGP>;
            const Eigen::Matrix<double, 1, NGP* NGP> e12_gauss = E12.row(i) * PSI<DGs, NGP>;
            const Eigen::Matrix<double, 1, NGP* NGP> e22_gauss = E22.row(i) * PSI<DGs, NGP>;

            Eigen::Matrix<double, 1, NGP* NGP> s11_gauss = S11.row(i) * PSI<DGs, NGP>;
            Eigen::Matrix<double, 1, NGP* NGP> s12_gauss = S12.row(i) * PSI<DGs, NGP>;
            Eigen::Matrix<double, 1, NGP* NGP> s22_gauss = S22.row(i) * PSI<DGs, NGP>;

            //! Current normal stress for the evaluation of tildeP (Eqn. 1)
            Eigen::Matrix<double, 1, NGP* NGP> sigma_n
                = 0.5 * (s11_gauss.array() + s22_gauss.array());

            //! exp(-C(1-A))
            const auto expC = (params.compaction_param * (1.0 - a_gauss.array())).exp();

            // Eqn. 25
            const auto powalphaexpC
                = ((1. - d_gauss.array()) * expC.array()).pow(params.exponent_relaxation_sigma - 1);
            const Eigen::Matrix<double, 1, NGP* NGP> time_viscous
                = params.undamaged_time_relaxation_sigma * powalphaexpC;

            //! BBM  Computing tildeP according to (Eqn. 7b and Eqn. 8)
            // (Eqn. 8)
            const Eigen::Matrix<double, 1, NGP* NGP> Pmax
                = params.P0 * h_gauss.array().pow(params.exponent_compression_factor) * expC;

            // (Eqn. 7b) Prepare tildeP
            // tildeP must be capped at 1 to get an elastic response
            // (Eqn. 7b) Select case based on sigma_n
            const Eigen::Matrix<double, 1, NGP* NGP> tildeP
                = (sigma_n.array() < 0.0)
                      .select((-Pmax.array() / sigma_n.array()).min(1.0).matrix(), 0.);

            // multiplicator
            const Eigen::Matrix<double, 1, NGP* NGP> multiplicator
                = time_viscous.array() / (time_viscous.array() + (1. - tildeP.array()) * dt_mom);

            //! Eqn. 9
            const Eigen::Matrix<double, 1, NGP* NGP> elasticity
                = (params.young * (1. - d_gauss.array()) * expC).matrix();

            // Eqn. 12: first factor on RHS
            /* Stiffness matrix
             * / (K:e)11 \       1     /  1  nu    0  \ / e11 \
             * | (K:e)22 |  =  ------- | nu   1    0  | | e22 |
             * \ (K:e)12 /    1 - nu^2 \  0   0  1-nu / \ e12 /
             */

            const Eigen::Matrix<double, 1, NGP* NGP> Dunit_factor
                = dt_mom * elasticity.array() / (1. - (params.nu0 * params.nu0));

            s11_gauss.array()
                += Dunit_factor.array() * (e11_gauss.array() + params.nu0 * e22_gauss.array());
            s22_gauss.array()
                += Dunit_factor.array() * (params.nu0 * e11_gauss.array() + e22_gauss.array());
            s12_gauss.array() += Dunit_factor.array() * e12_gauss.array() * (1. - params.nu0);

            //! Implicit part of RHS (Eqn. 33)
            s11_gauss.array() *= multiplicator.array();
            s22_gauss.array() *= multiplicator.array();
            s12_gauss.array() *= multiplicator.array();

            sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());
            const Eigen::Matrix<double, 1, NGP* NGP> tau
                = (0.25 * (s11_gauss.array() - s22_gauss.array()).square()
                    + s12_gauss.array().square())
                      .sqrt();

            const double scale_coef = std::sqrt(0.1 / smesh.h(i));
            const double cohesion = params.C_lab * scale_coef;
            const double compr_strength = params.compr_strength * scale_coef;

            // Mohr-Coulomb failure using Mssrs. Plante & Tremblay's formulation
            // sigma_s + tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
            Eigen::Matrix<double, 1, NGP* NGP> dcrit
                = (tau.array() + params.tan_phi * sigma_n.array() > 0.)
                      .select(cohesion / (tau.array() + params.tan_phi * sigma_n.array()), 1.);

            // Compressive failure using Mssrs. Plante & Tremblay's formulation
            dcrit = (sigma_n.array() < -compr_strength)
                        .select(-compr_strength / sigma_n.array(), dcrit);

            // Only damage when we're outside
            dcrit = dcrit.array().min(1.0);

            const Eigen::Matrix<double, 1, NGP* NGP> td = smesh.h(i)
                * std::sqrt(2. * (1. + params.nu0) * params.rho_ice) / elasticity.array().sqrt();

            // Update damage
            d_gauss.array() += (1. - d_gauss.array()) * (1. - dcrit.array()) * dt_mom / td.array();

            // Relax stress in Gassus points
            s11_gauss.array() -= s11_gauss.array() * (1. - dcrit.array()) * dt_mom / td.array();
            s12_gauss.array() -= s12_gauss.array() * (1. - dcrit.array()) * dt_mom / td.array();
            s22_gauss.array() -= s22_gauss.array() * (1. - dcrit.array()) * dt_mom / td.array();

            // INTEGRATION OF STRESS AND DAMAGE
            const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> J
                = ParametricTools::J<3>(smesh, i);
            // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // with the gauss weights and with J. This is a 8 x 9 matrix
            const Eigen::Matrix<Nextsim::FloatType, DGs, NGP* NGP> imass_psi = ParametricTools::massMatrix<DGs>(smesh, i).inverse()
                * (PSI<DGs, NGP>.array().rowwise() * (GAUSSWEIGHTS<NGP>.array() * J.array())).matrix();

            S11.row(i) = imass_psi * s11_gauss.matrix().transpose();
            S12.row(i) = imass_psi * s12_gauss.matrix().transpose();
            S22.row(i) = imass_psi * s22_gauss.matrix().transpose();

            const Eigen::Matrix<Nextsim::FloatType, DGa, NGP* NGP> imass_psi2 = ParametricTools::massMatrix<DGa>(smesh, i).inverse()
                * (PSI<DGa, NGP>.array().rowwise() * (GAUSSWEIGHTS<NGP>.array() * J.array())).matrix();

            D.row(i) = imass_psi2 * d_gauss.matrix().transpose();
        }
    }

#undef NGP

} /* namespace MEB */

} /* namespace Nextsim */

#endif /* __MEB_HPP */
