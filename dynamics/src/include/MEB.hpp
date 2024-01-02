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



namespace Nextsim {

/*!
 * This namespace collects the routines required for the MEB solver
 */
namespace MEB {

    inline constexpr double SQR(double x) { return x * x; }

    /*!
     * @brief Calculate Stresses for the current time step and update damage.
     *
     * @details MEB model, numbering of equations in comments according to:
     *          Dansereau et al. 2016
     *          "A Maxwell elasto-brittle rheology for sea ice modelling"
     *          https://doi:10.5194/tc-10-1339-2016
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
        constexpr size_t NGaussPoints = (DGs == 8 ? 3 : (DGs == 3 ? 2 : -1));
        typedef Eigen::Matrix<double, 1, NGaussPoints * NGaussPoints> GaussPointVector;


        //! Stress and Damage Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            //! Evaluate values in Gauss points (3 point Gauss rule in 2d => 9 points)
            const GaussPointVector h_gauss = (H.row(i) * PSI<DGa, NGaussPoints>).array().max(0.0).matrix();
            const GaussPointVector a_gauss = (A.row(i) * PSI<DGa, NGaussPoints>).array().max(0.0).min(1.0).matrix();
            GaussPointVector d_gauss = (D.row(i) * PSI<DGa, NGaussPoints>).array().max(1e-12).min(1.0).matrix();

            const GaussPointVector e11_gauss = E11.row(i) * PSI<DGs, NGaussPoints>;
            const GaussPointVector e12_gauss = E12.row(i) * PSI<DGs, NGaussPoints>;
            const GaussPointVector e22_gauss = E22.row(i) * PSI<DGs, NGaussPoints>;

            GaussPointVector s11_gauss = S11.row(i) * PSI<DGs, NGaussPoints>;
            GaussPointVector s12_gauss = S12.row(i) * PSI<DGs, NGaussPoints>;
            GaussPointVector s22_gauss = S22.row(i) * PSI<DGs, NGaussPoints>;

            //! exp(-C(1-A))
            const GaussPointVector expC = (params.compaction_param * (1.0 - a_gauss.array())).exp().array();

            // Eqn. 20
            GaussPointVector powalpha = (d_gauss.array()).pow(params.exponent_relaxation_sigma - 1.);
            const GaussPointVector time_viscous = (params.undamaged_time_relaxation_sigma * powalpha.array() ).matrix();

            // Eqn. 4: first factor on RHS
            const double Dunit_factor = 1. / (1. - (params.nu0 * params.nu0));

            //! MEB
            // 1. / (1. + dt / lambda) 
            GaussPointVector multiplicator = (1. / (1. + dt_mom / time_viscous.array())).matrix();

            
            //! Eqn. 24
            const GaussPointVector elasticity = h_gauss.array() * params.young * d_gauss.array() * expC.array();

            // Eqn. 4: first factor on RHS
            /* Stiffness matrix
             * / (K:e)11 \       1     /  1  nu    0  \ / e11 \
             * | (K:e)22 |  =  ------- | nu   1    0  | | e22 |
             * \ (K:e)12 /    1 - nu^2 \  0   0  1-nu / \ e12 /
             */


            s11_gauss += (dt_mom * 1. / (1. + params.nu0) * (elasticity.array() * e11_gauss.array())).matrix()
                + (dt_mom * Dunit_factor * params.nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array()))).matrix();
            s12_gauss += (dt_mom * 1. / (1. + params.nu0) * (elasticity.array() * e12_gauss.array())).matrix();
            s22_gauss += (dt_mom * 1. / (1. + params.nu0) * (elasticity.array() * e22_gauss.array())).matrix()
                + (dt_mom * Dunit_factor * params.nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array()))).matrix();


            //! Implicit part of RHS (Eqn. 3)
            s11_gauss.array() *= multiplicator.array();
            s12_gauss.array() *= multiplicator.array();
            s22_gauss.array() *= multiplicator.array();

            //! Current normal and tangent stress for the evaluation Mohr-Couloumb (Eqn 13)
            const GaussPointVector sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());
            const GaussPointVector tau = (0.25 * (s11_gauss.array() - s22_gauss.array()).square() + s12_gauss.array().square()).sqrt();

            GaussPointVector dcrit = GaussPointVector::Ones();

            // Fixed Cohesion
            const GaussPointVector cohesion = params.c0 * h_gauss.array() ;
            
            
            //! This is not part of Dansereau et al. 2016
            const double scale_coef = std::sqrt(0.1 / smesh.h(i));
            const double compr_strength = params.compr_strength * scale_coef;

            // Mohr-Coulomb failure using Mssrs. Plante & Tremblay's formulation
            // sigma_s + tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
            dcrit
                = (tau.array() + params.tan_phi * sigma_n.array() > 0.)
                      .select(cohesion.array() / (tau.array() + params.tan_phi * sigma_n.array()), 1.);

            // Compressive failure using Mssrs. Plante & Tremblay's formulation
            dcrit = (sigma_n.array() < -compr_strength)
                        .select(-compr_strength / sigma_n.array(), dcrit);
            // Only damage when we're outside
            dcrit = dcrit.array().min(1.0);

            const double C_e = 500.0; //!  Elatic shear wave propagation speed in m/s, see Table 1. and Section 4.1.1 
            const double td = smesh.h(i)/C_e  ;

            // Update damage
            d_gauss.array() -= d_gauss.array() * (1. - dcrit.array()) * dt_mom / td;

            // Relax stress in Gassus points
            s11_gauss.array() -= s11_gauss.array() * (1. - dcrit.array()) * dt_mom / td;
            s12_gauss.array() -= s12_gauss.array() * (1. - dcrit.array()) * dt_mom / td;
            s22_gauss.array() -= s22_gauss.array() * (1. - dcrit.array()) * dt_mom / td;


            // INTEGRATION OF STRESS AND DAMAGE
            const Eigen::Matrix<Nextsim::FloatType, 1, NGaussPoints * NGaussPoints> J = ParametricTools::J<3>(smesh, i);
            // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // with the gauss weights and with J. This is a 8 x 9 matrix
            const Eigen::Matrix<Nextsim::FloatType, DGs, NGaussPoints * NGaussPoints> imass_psi = ParametricTools::massMatrix<DGs>(smesh, i).inverse()
                * (PSI<DGs, NGaussPoints>.array().rowwise() * (GAUSSWEIGHTS<NGaussPoints>.array() * J.array())).matrix();

            S11.row(i) = imass_psi * s11_gauss.matrix().transpose();
            S12.row(i) = imass_psi * s12_gauss.matrix().transpose();
            S22.row(i) = imass_psi * s22_gauss.matrix().transpose();

            const Eigen::Matrix<Nextsim::FloatType, DGa, NGaussPoints * NGaussPoints> imass_psi2 = ParametricTools::massMatrix<DGa>(smesh, i).inverse()
                * (PSI<DGa, NGaussPoints>.array().rowwise() * (GAUSSWEIGHTS<NGaussPoints>.array() * J.array())).matrix();

            D.row(i) = imass_psi2 * d_gauss.matrix().transpose();
        }
    }

} /* namespace MEB */

} /* namespace Nextsim */

#endif /* __MEB_HPP */
