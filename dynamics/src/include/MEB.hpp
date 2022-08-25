/*!
 * @file MEB.hpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MEB_HPP
#define __MEB_HPP

//#include <Eigen/Dense>

#include "MEBParameters.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace RefScaleMain {

constexpr double rho_ice = 900.0; //!< Sea ice density
constexpr double rho_atm = 1.3; //!< Air density
constexpr double rho_ocean = 1026.0; //!< Ocean density

constexpr double C_atm = 1.2e-3; //!< Air drag coefficient
constexpr double C_ocean = 5.5e-3; //!< Ocean drag coefficient

constexpr double F_atm = C_atm * rho_atm; //!< effective factor for atm-forcing
constexpr double F_ocean = C_ocean * rho_ocean; //!< effective factor for ocean-forcing

// parameters form nextsim options.cpp line 302
constexpr double compaction_param = -20.; //!< Compation parameter

constexpr double nu0 = 1. / 3.; //!< \param Poisson's ratio
}

// Benchmark testcase from [Plante Tremblay, 2021]
namespace RefScaleCanada {
using namespace RefScaleMain;

constexpr double T = 4 * 60. * 60.; //!< Time horizon 2 hours
constexpr double L = 60000.0; //!< Size of domain
constexpr double vmax_ocean = 0.0; //!< Maximum velocity of ocean
constexpr double vmax_atm = 20; //!< Max. vel. of wind

// Parameters from Table 1. [Plante Tremblay, 2021]
constexpr double damage_timescale = 1.; //<! Damage timescale
constexpr double undamaged_time_relaxation_sigma = 1e5; //!< Test more viscous
constexpr double exponent_relaxation_sigma = 3;
constexpr double young = 1e9;
constexpr double compr_strength
    = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]

// TODO missing 45\deg it goes to Compresssion
constexpr double tan_phi = 0.7; //! \param tan_phi (double) Internal friction coefficient (mu)
constexpr double sin_phi = 0.7; //! \param sin_phi (double) Internal friction coefficient (mu)

constexpr double c0 = 10.e3; //! \param
constexpr double sigma_c0 = 50.e3; //! \param

// constexpr double C_lab = 2.0e6; //! \param C_lab (double) Test [Pa]
// constexpr double time_relaxation_damage = 2160000.; //!< 25 days in seconds
//
// constexpr double compression_factor = 10e3; //! \param Max pressure for damaged converging ice
// constexpr double exponent_compression_factor
//     = 1.5; //! \param Power of ice thickness in the pressure coefficient

constexpr double fc = 0.0; //!< Coriolis

}

// Benchmark testcase from [Mehlmann / Richter, ...]
namespace RefScale {
using namespace RefScaleMain;

constexpr double T = 2 * 24. * 60. * 60.; //!< Time horizon 2 days
constexpr double L = 512000.0; //!< Size of domain
constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
constexpr double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind

constexpr double Pstar = 27500.0; //!< Ice strength
constexpr double fc = 1.46e-4; //!< Coriolis
constexpr double DeltaMin = 2.e-9; //!< Viscous regime

// parameters form nextsim options.cpp line 302
constexpr double compaction_param = -20.; //!< Compation parameter

// constexpr double undamaged_time_relaxation_sigma = 1e7; //!< seconds
constexpr double undamaged_time_relaxation_sigma = 1e5; //!< Test more viscous

constexpr double time_relaxation_damage = 2160000.; //!< 25 days in seconds
constexpr double compression_factor = 10e3; //! \param Max pressure for damaged converging ice
constexpr double exponent_compression_factor
    = 1.5; //! \param Power of ice thickness in the pressure coefficient

constexpr double exponent_relaxation_sigma = 5;
constexpr double young = 5.9605e+08;
constexpr double nu0 = 1. / 3.; //!< \param Poisson's ratio
constexpr double compr_strength
    = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]
constexpr double tan_phi = 0.7; //! \param tan_phi (double) Internal friction coefficient (mu)
// constexpr double C_lab = 2.0e6; //! \param C_lab (double) Cohesion at the lab scale (10 cm) [Pa]
constexpr double C_lab = 2.0e6; //! \param C_lab (double) Test [Pa]

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
     * @details BBM model
     *
     * @tparam DGstress Stress adn Strain DG degree
     * @tparam DGtracer H, A, D DG degree
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
     * @param dt_momentum timestep for momentum subcycle
     */
    template <int DGstress, int DGtracer>
    void StressUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, CellVector<DGtracer>& D, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // There's no ice so we set sigma to 0 and carry on
            const double min_c = 0.1;
            if (A(i, 0) <= min_c) {
                D(i, 0) = 0.0;
                S11.row(i) *= 0.0;
                S12.row(i) *= 0.0;
                S22.row(i) *= 0.0;
                continue;
            }

            //! - Updates the internal stress
            double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
            double const time_viscous = RefScale::undamaged_time_relaxation_sigma
                * std::pow((1. - D(i, 0)) * expC, RefScale::exponent_relaxation_sigma - 1.);

            //! Plastic failure tildeP
            double tildeP;
            if (sigma_n < 0.) {
                double const Pmax = RefScale::compression_factor
                    * pow(H(i, 0), RefScale::exponent_compression_factor)
                    * exp(RefScale::compaction_param * (1.0 - A(i, 0)));
                // tildeP must be capped at 1 to get an elastic response
                tildeP = std::min(1., -Pmax / sigma_n);
            } else {
                tildeP = 0.;
            }

            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
            // double const multiplicator
            //    = std::min(1. - 1e-12, time_viscous / (time_viscous + dt_momentum * (1. - tildeP)));

            double const multiplicator = 1. / (1. + dt_momentum / time_viscous * (1. - tildeP));
            tildeP = multiplicator;

            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;
            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * elasticity
                * (1 / (1 + RefScale::nu0) * E11.row(i)
                    + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += dt_momentum * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
            S22.row(i) += dt_momentum * elasticity
                * (1 / (1 + RefScale::nu0) * E22.row(i)
                    + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S11.row(i) *= multiplicator;
            S12.row(i) *= multiplicator;
            S22.row(i) *= multiplicator;

            //! - Estimates the level of damage from updated stress and the local damage criterion

            // continiue if stress in inside the failure envelope
            //  Compute the shear and normal stresses, which are two invariants of the internal
            //  stress tensor
            double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
            // update sigma_n
            sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));

            // Cohesion Eqn. (21)
            // Reference length scale is fixed 0.1 since its cohesion parameter at the lab scale (10
            // cm)
            double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);

            // d critical Eqn. (29)
            double dcrit;
            if (sigma_n < -RefScale::compr_strength)
                dcrit = -RefScale::compr_strength / sigma_n;
            else
                // M_Cohesion[cpt] depends on local random contribution
                // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]);
                // finiteelement.cpp#L3834
                dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            // Calcutate d_critical as Veronique:

            // Calculate the characteristic time for damage and damage increment
            // M_delta_x[cpt] = mesh.hx ???
            double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                / std::sqrt(elasticity);
            // Calculate the adjusted level of damage
            if ((0. < dcrit) && (dcrit < 1.)) {
                // Eqn. (34)
                D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                // Recalculate the new state of stress by relaxing elstically Eqn. (36)
                S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
            }

            // Relax damage
            D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage * std::exp(RefScale::compaction_param * (1. - A(i, 0))));
            // D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage * std::exp(RefScale::compaction_param * (1. - A(i, 0))));
        }
    }

    /*ELO /

        //! DG0 for Stress and Strain
        template <int DGstress, int DGtracer>
        void StressUpdateMEB(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
            CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
            const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
            const CellVector<DGtracer>& A, CellVector<DGtracer>& D, const double dt_momentum)
    {
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            //======================================================================
            //! Visco-Elasatic Prediction
            //======================================================================
            double const expC = std::exp(RefScaleCanada::compaction_param * (1. - A(i, 0)));
            double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));

            // Eqn. 25
            // double time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma
            //    * std::pow((1. - D(i, 0)), RefScaleCanada::exponent_relaxation_sigma - 1.) / (H(i, 0) * expC);
            // Eqn 9.
            double time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma
                * std::pow((1. - D(i, 0)), RefScaleCanada::exponent_relaxation_sigma - 1.);

            // Eqn. 24 Additional multiplic0cation by H
            double elasticity = RefScaleCanada::young * H(i, 0) * (1. - D(i, 0)) * expC;

            // 1. / (1. + dt / lambda) Eqn. 18
            double multiplicator = 1. / (1. + dt_momentum / time_viscous);

            time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma;
            // elasticity = RefScaleCanada::young;

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScaleCanada::nu0) * E11.row(i)
                    + Dunit_factor * RefScaleCanada::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += dt_momentum * elasticity * 1. / (1. + RefScaleCanada::nu0) * E12.row(i);
            S22.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScaleCanada::nu0) * E22.row(i)
                    + Dunit_factor * RefScaleCanada::nu0 * (E11.row(i) + E22.row(i)));

            // BBM
            double const Pmax = RefScale::compression_factor * pow(H(i, 0), 1.5); // * exp(-20.0 * (1.0 - A(i, 0)));
            // Strange Discountinious function
            if (sigma_n < 0.) {
                double tildeP = std::min(1., -Pmax / sigma_n);
                // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
                multiplicator = 1. / (1. + (1. - tildeP) * dt_momentum / time_viscous);
            }

            S12.row(i) *= multiplicator;
            S11.row(i) *= multiplicator;
            S22.row(i) *= multiplicator;

            //======================================================================
            //! - Estimates the level of damage from the updated internal stress and the local
            //! damage criterion
            //======================================================================

            // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
            sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            // shear stress
            double tau = std::sqrt(0.25 * (S11(i, 0) - S22(i, 0)) * (S11(i, 0) - S22(i, 0)) + S12(i, 0) * S12(i, 0));
            assert(tau >= 0);

            //
            const double c = RefScaleCanada::c0 * H(i, 0) * expC;
            const double sigma_c = RefScaleCanada::sigma_c0 * H(i, 0) * expC;

            double dcrit(1.0);
            // Regime(i, 0) = 0;
            //// Mohr-Coulomb condition
            // if (tau + RefScaleCanada::sin_phi * sigma_n - c < 0) {
            //     Regime(i, 0) = 1;
            // }
            // if (sigma_n - tau > sigma_c) {
            //     Regime(i, 0) = 2;
            // }

            if (tau + RefScaleCanada::sin_phi * sigma_n - c > 0)
                dcrit = std::min(1., c / (tau + RefScaleCanada::sin_phi * sigma_n));

            // Relax stress
            S11.row(i) *= dcrit;
            S12.row(i) *= dcrit;
            S22.row(i) *= dcrit;

            // Update damage
            D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / RefScaleCanada::damage_timescale;
        }
    }

    //! CG1-DG1
    void StressUpdateMEB(const Mesh& mesh, CellVector<3>& S11, CellVector<3>& S12,
        CellVector<3>& S22, const CellVector<3>& E11, const CellVector<3>& E12,
        const CellVector<3>& E22, const CellVector<3>& H,
        const CellVector<3>& A, CellVector<3>& D, const double dt_momentum, CellVector<3>& Dcrit)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            //!  number of gauss points in each direction equalls 3, for 2d 3^2
            const Eigen::Matrix<double, 1, 9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            // const Eigen::Matrix<double, 1, 9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();
            const Eigen::Matrix<double, 1, 9> d_gauss = (D.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const Eigen::Matrix<double, 1, 9> e11_gauss = E11.block<1, 3>(i, 0) * PSI<3,3>;
            const Eigen::Matrix<double, 1, 9> e12_gauss = E12.block<1, 3>(i, 0) * PSI<3,3>;
            const Eigen::Matrix<double, 1, 9> e22_gauss = E22.block<1, 3>(i, 0) * PSI<3,3>;

            Eigen::Matrix<double, 1, 9> s11_gauss = S11.block<1, 3>(i, 0) * PSI<3,3>;
            Eigen::Matrix<double, 1, 9> s12_gauss = S12.block<1, 3>(i, 0) * PSI<3,3>;
            Eigen::Matrix<double, 1, 9> s22_gauss = S22.block<1, 3>(i, 0) * PSI<3,3>;

            auto expC = 1.; //(-20.0 * (1.0 - a_gauss.array())).exp();
            // Eqn. 24 Additional multiplic0cation by H
            // const Eigen::Matrix<double, 1, 9> elasticity = (RefScaleCanada::young * h_gauss.array() * (1. - d_gauss.array()) * expC.array()).matrix();
            const Eigen::Matrix<double, 1, 9> elasticity = (RefScaleCanada::young * h_gauss.array() * (1. - d_gauss.array())).matrix();

            auto powalpha = (1. - d_gauss.array()).pow(RefScaleCanada::exponent_relaxation_sigma - 1.);
            // Eqn. 25
            const Eigen::Matrix<double, 1, 9> time_viscous = (RefScaleCanada::undamaged_time_relaxation_sigma * powalpha.array()).matrix();

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // 1. / (1. + dt / lambda) Eqn. 18
            Eigen::Matrix<double, 1, 9> multiplicator = (1. / (1. + dt_momentum / time_viscous.array())).matrix();

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e11_gauss.array()).matrix() * IBC33
                + dt_momentum * Dunit_factor * RefScaleCanada::nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array())).matrix() * IBC33;

            S12.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e12_gauss.array()).matrix() * IBC33;

            S22.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e22_gauss.array()).matrix() * IBC33
                + dt_momentum * Dunit_factor * RefScaleCanada::nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array())).matrix() * IBC33;

            // BBM
            Eigen::Matrix<double, 1, 9> sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());

            const Eigen::Matrix<double, 1, 9> Pmax = RefScale::compression_factor * h_gauss.array().pow(1.5);
            // Strange Discountinious function
            const Eigen::Matrix<double, 1, 9> tildeP = (-Pmax.array() / sigma_n.array()).min(1.0).matrix();

            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
            multiplicator = (sigma_n.array() < 0.).select((1. / (1. + (1. - tildeP.array()) * dt_momentum / time_viscous.array())).matrix(), multiplicator);

            // multiplicator = multiplicator_gp * IBC33;
            //  Eqn. (34)
            //! Coefficient-wise multiplication
            S11.row(i) = (s11_gauss.array() * multiplicator.array()).matrix() * IBC33;
            S12.row(i) = (s12_gauss.array() * multiplicator.array()).matrix() * IBC33;
            S22.row(i) = (s22_gauss.array() * multiplicator.array()).matrix() * IBC33;

            //======================================================================
            //! - Estimates the level of damage from the updated internal stress and the local
            //! damage criterion
            //======================================================================

            s11_gauss = S11.block<1, 3>(i, 0) * PSI<3,3>;
            s12_gauss = S12.block<1, 3>(i, 0) * PSI<3,3>;
            s22_gauss = S22.block<1, 3>(i, 0) * PSI<3,3>;

            // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
            // double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());

            // shear stress
            // double tau = std::sqrt(0.25 * (S11(i, 0) - S22(i, 0)) * (S11(i, 0) - S22(i, 0)) + S12(i, 0) * S12(i, 0));

            const Eigen::Matrix<double, 1, 9> tau = (0.25 * (s11_gauss.array() - s22_gauss.array()).square() + s12_gauss.array().square()).sqrt();
            // assert(tau >= 0);

            // const double c = RefScaleCanada::c0 * H(i, 0) * std::exp(RefScaleCanada::compaction_param * (1. - A(i, 0)));
            // const double sigma_c = RefScaleCanada::sigma_c0 * H(i, 0) * std::exp(RefScaleCanada::compaction_param * (1. - A(i, 0)));
            const Eigen::Matrix<double, 1, 9> c = RefScaleCanada::c0 * h_gauss.array() * expC;
            // const Eigen::Matrix<double, 1, 9> sigma_c = RefScaleCanada::sigma_c0 * h_gauss.array() * expC;

            // double dcrit(1.0);
            Eigen::Matrix<double, 1, 9> dcrit;
            dcrit << 1., 1., 1., 1., 1., 1., 1., 1., 1.;
            // if (tau + RefScaleCanada::sin_phi * sigma_n - c > 0)
            //     dcrit = std::min(1., c / (tau + RefScaleCanada::sin_phi * sigma_n));

            // Mohr-Coulomb Criterion
            dcrit = (tau.array() + RefScaleCanada::sin_phi * sigma_n.array() - c.array() > 0).select(c.array() / (tau.array() + RefScaleCanada::sin_phi * sigma_n.array()), dcrit);
            // Compression cutoff
            // dcrit = (sigma_c.array() / (sigma_n.array() - tau.array()) < dcrit.array()).select(sigma_c.array() / (sigma_n.array() - tau.array()), dcrit);

            dcrit = dcrit.array().min(1.0).matrix();

            auto dcrit_int = dcrit * IBC33;
            Dcrit.row(i) = dcrit_int;

            // S11.row(i).array() *= dcrit_int.array();
            // S12.row(i).array() *= dcrit_int.array();
            // S22.row(i).array() *= dcrit_int.array();

            // Relax stress in Gassus points
            S11.row(i) = (s11_gauss.array() * dcrit.array()).matrix() * IBC33;
            S12.row(i) = (s12_gauss.array() * dcrit.array()).matrix() * IBC33;
            S22.row(i) = (s22_gauss.array() * dcrit.array()).matrix() * IBC33;

            // Update damage
            // D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / RefScaleCanada::damage_timescale;
            D.row(i) += ((1.0 - d_gauss.array()) * (1.0 - dcrit.array()) * dt_momentum / RefScaleCanada::damage_timescale).matrix() * IBC33;
        }
    }
    ELO */

    /*!
     * @brief Calculate Stresses for the current time step and update damage.
     *
     * @details BBM model, numbering of equations in comments according to:
     *          https://doi.org/10.1029/2021MS002685
     *
     * @tparam DGstress Stress adn Strain DG degree
     * @tparam DGtracer H, A, D DG degree
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
    void StressUpdateHighOrder(const MEBParameters& params,
        const SasipMesh& smesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<3>& H,
        const CellVector<3>& A, CellVector<3>& D,
        const double dt_mom)
    {

        //! Stress and Damage Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            //! Evaluate values in Gauss points (3 point Gauss rule in 2d => 9 points)
            const Eigen::Matrix<double, 1, 9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3, 3>).array().max(0.0).matrix();
            const Eigen::Matrix<double, 1, 9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3, 3>).array().max(0.0).min(1.0).matrix();
            Eigen::Matrix<double, 1, 9> d_gauss = (D.block<1, 3>(i, 0) * PSI<3, 3>).array().max(0.0).min(1.0).matrix();

            const Eigen::Matrix<double, 1, 9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8, 3>;
            const Eigen::Matrix<double, 1, 9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8, 3>;
            const Eigen::Matrix<double, 1, 9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8, 3>;

            Eigen::Matrix<double, 1, 9> s11_gauss = S11.block<1, 8>(i, 0) * PSI<8, 3>;
            Eigen::Matrix<double, 1, 9> s12_gauss = S12.block<1, 8>(i, 0) * PSI<8, 3>;
            Eigen::Matrix<double, 1, 9> s22_gauss = S22.block<1, 8>(i, 0) * PSI<8, 3>;

            //! Current normal stress for the evaluation of tildeP (Eqn. 1)
            Eigen::Matrix<double, 1, 9> sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());

            auto expC = (-20.0 * (1.0 - a_gauss.array())).exp();

            const Eigen::Matrix<double, 1, 9> elasticity = (params.young * h_gauss.array() * (1. - d_gauss.array())).matrix();

            auto powalpha = (1. - d_gauss.array()).pow(params.exponent_relaxation_sigma - 1.);

            // Eqn. 25
            const Eigen::Matrix<double, 1, 9> time_viscous = (params.undamaged_time_relaxation_sigma * powalpha.array()).matrix();

            // Eqn. 12: first factor on RHS
            double const Dunit_factor = 1. / (1. - (params.nu0 * params.nu0));

            //! MEB
            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 34
            // 1. / (1. + dt / lambda) Eqn. 33-34
            Eigen::Matrix<double, 1, 9> multiplicator = (1. / (1. + dt_mom / time_viscous.array())).matrix();

            //! BBM  Computing tildeP according to (Eqn. 7b and Eqn. 8)
            // (Eqn. 8)
            const Eigen::Matrix<double, 1, 9> Pmax = params.P0 * h_gauss.array().pow(1.5);

            // (Eqn. 7b) Prepare tildeP
            const Eigen::Matrix<double, 1, 9> tildeP = (-Pmax.array() / sigma_n.array()).min(1.0).matrix();
            // (Eqn. 7b) Select case based on sigma_n
            // multiplicator = (sigma_n.array() < 0.).select((1. / (1. + (1. - tildeP.array()) * dt_mom / time_viscous.array())).matrix(), multiplicator);

            //! Implicit part of RHS (Eqn. 33)
            s11_gauss.array() *= multiplicator.array();
            s12_gauss.array() *= multiplicator.array();
            s22_gauss.array() *= multiplicator.array();

            sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());
            const Eigen::Matrix<double, 1, 9> tau = (0.25 * (s11_gauss.array() - s22_gauss.array()).square() + s12_gauss.array().square()).sqrt();

            const Eigen::Matrix<double, 1, 9> c = params.c0 * h_gauss.array() * expC;

            Eigen::Matrix<double, 1, 9> dcrit = { 1., 1., 1., 1., 1., 1., 1., 1., 1. };

            // Mohr-Coulomb Criterion
            dcrit = (tau.array() + params.sin_phi * sigma_n.array() - c.array() > 0).select(c.array() / (tau.array() + params.sin_phi * sigma_n.array()), dcrit);
            // Compression cutoff
            // dcrit = (sigma_c.array() / (sigma_n.array() - tau.array()) < dcrit.array()).select(sigma_c.array() / (sigma_n.array() - tau.array()), dcrit);

            dcrit = dcrit.array().min(1.0);

            // Relax stress in Gassus points
            s11_gauss.array() *= dcrit.array();
            s12_gauss.array() *= dcrit.array();
            s22_gauss.array() *= dcrit.array();

            // Update damage
            // D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / RefScaleCanada::damage_timescale;
            d_gauss = ((1.0 - d_gauss.array()) * (1.0 - dcrit.array()) * dt_mom / params.damage_timescale);

            // INTEGRATION OF STRESS AND DAMAGE
            const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
            // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // with the gauss weights and with J. This is a 8 x 9 matrix
            const Eigen::Matrix<Nextsim::FloatType, 8, 9> imass_psi = ParametricTools::massMatrix<8>(smesh, i).inverse()
                * (PSI<8, 3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J.array())).matrix();

            S11.row(i) = imass_psi * s11_gauss.matrix().transpose();
            S12.row(i) = imass_psi * s12_gauss.matrix().transpose();
            S22.row(i) = imass_psi * s22_gauss.matrix().transpose();

            const Eigen::Matrix<Nextsim::FloatType, 3, 9> imass_psi2 = ParametricTools::massMatrix<3>(smesh, i).inverse()
                * (PSI<3, 3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J.array())).matrix();

            D.row(i) += imass_psi2 * d_gauss.matrix().transpose();
        }
    }

} /* namespace MEB */

} /* namespace Nextsim */

#endif /* __MEB_HPP */
