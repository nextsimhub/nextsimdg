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

    void StressUpdateHighOrder(const MEBParameters& vpparameters,
        const SasipMesh& smesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<3>& H,
        const CellVector<3>& A, CellVector<3>& D,
        const double dt_mom)
    {

        double alpha = 800;
        double beta = 800;

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const Eigen::Matrix<double, 1, 9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const Eigen::Matrix<double, 1, 9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();
            const Eigen::Matrix<double, 1, 9> d_gauss = (D.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const Eigen::Matrix<double, 1, 9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8,3>;
            const Eigen::Matrix<double, 1, 9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8,3>;
            const Eigen::Matrix<double, 1, 9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8,3>;

            Eigen::Matrix<double, 1, 9> s11_gauss = S11.block<1, 8>(i, 0) * PSI<8,3>;
            Eigen::Matrix<double, 1, 9> s12_gauss = S12.block<1, 8>(i, 0) * PSI<8,3>;
            Eigen::Matrix<double, 1, 9> s22_gauss = S22.block<1, 8>(i, 0) * PSI<8,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

            // std::cout << P << std::endl << DELTA << std::endl;
            // abort();

            auto expC = 1.; //(-20.0 * (1.0 - a_gauss.array())).exp();

            const Eigen::Matrix<double, 1, 9> elasticity = (vpparameters.young * h_gauss.array() * (1. - d_gauss.array())).matrix();

            auto powalpha = (1. - d_gauss.array()).pow(vpparameters.exponent_relaxation_sigma - 1.);
            // Eqn. 25
            const Eigen::Matrix<double, 1, 9> time_viscous = (vpparameters.undamaged_time_relaxation_sigma * powalpha.array()).matrix();

            double const Dunit_factor = 1. / (1. - (vpparameters.nu0 * vpparameters.nu0));

            // 1. / (1. + dt / lambda) Eqn. 18
            Eigen::Matrix<double, 1, 9> multiplicator = (1. / (1. + dt_mom / time_viscous.array())).matrix();

            // Elasit prediction Eqn. (32)
            const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
            // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // with the gauss weights and with J. This is a 8 x 9 matrix
            const Eigen::Matrix<Nextsim::FloatType, 8, 9> imass_psi = ParametricTools::massMatrix<8>(smesh, i).inverse()
                * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS_3.array() * J.array())).matrix();

            S11.row(i) += imass_psi * (dt_mom * 1. / (1. + vpparameters.nu0) * (elasticity.array() * e11_gauss.array())).matrix().transpose()
                + imass_psi * (dt_mom * Dunit_factor * vpparameters.nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array()))).matrix().transpose();

            S12.row(i) += imass_psi * (dt_mom * 1. / (1. + vpparameters.nu0) * (elasticity.array() * e12_gauss.array())).matrix().transpose();

            S22.row(i) += imass_psi * (dt_mom * 1. / (1. + vpparameters.nu0) * (elasticity.array() * e22_gauss.array())).matrix().transpose()
                + imass_psi * (dt_mom * Dunit_factor * vpparameters.nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array()))).matrix().transpose();

            // BBM
            Eigen::Matrix<double, 1, 9> sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());

            const Eigen::Matrix<double, 1, 9> Pmax = RefScale::compression_factor * h_gauss.array().pow(1.5);
            // Strange Discountinious function
            const Eigen::Matrix<double, 1, 9> tildeP = (-Pmax.array() / sigma_n.array()).min(1.0).matrix();

            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
            multiplicator = (sigma_n.array() < 0.).select((1. / (1. + (1. - tildeP.array()) * dt_mom / time_viscous.array())).matrix(), multiplicator);

            // multiplicator = multiplicator_gp * IBC33;
            //  Eqn. (34)
            //! Coefficient-wise multiplication
            S11.row(i) = imass_psi * (s11_gauss.array() * multiplicator.array()).matrix().transpose();
            S12.row(i) = imass_psi * (s12_gauss.array() * multiplicator.array()).matrix().transpose();
            S22.row(i) = imass_psi * (s22_gauss.array() * multiplicator.array()).matrix().transpose();

            //======================================================================
            //! - Estimates the level of damage from the updated internal stress and the local
            //! damage criterion
            //======================================================================
            s11_gauss = S11.block<1, 8>(i, 0) * PSI<8,3>;
            s12_gauss = S12.block<1, 8>(i, 0) * PSI<8,3>;
            s22_gauss = S22.block<1, 8>(i, 0) * PSI<8,3>;

            sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());
            const Eigen::Matrix<double, 1, 9> tau = (0.25 * (s11_gauss.array() - s22_gauss.array()).square() + s12_gauss.array().square()).sqrt();

            const Eigen::Matrix<double, 1, 9> c = RefScaleCanada::c0 * h_gauss.array() * expC;

            Eigen::Matrix<double, 1, 9> dcrit;
            dcrit << 1., 1., 1., 1., 1., 1., 1., 1., 1.;

            // Mohr-Coulomb Criterion
            dcrit = (tau.array() + RefScaleCanada::sin_phi * sigma_n.array() - c.array() > 0).select(c.array() / (tau.array() + RefScaleCanada::sin_phi * sigma_n.array()), dcrit);
            // Compression cutoff
            // dcrit = (sigma_c.array() / (sigma_n.array() - tau.array()) < dcrit.array()).select(sigma_c.array() / (sigma_n.array() - tau.array()), dcrit);

            dcrit = dcrit.array().min(1.0).matrix();

            // auto dcrit_int = dcrit * IBC33;
            // Dcrit.row(i) = dcrit_int;

            // S11.row(i).array() *= dcrit_int.array();
            // S12.row(i).array() *= dcrit_int.array();
            // S22.row(i).array() *= dcrit_int.array();

            // Relax stress in Gassus points
            S11.row(i) = imass_psi * (s11_gauss.array() * dcrit.array()).matrix().transpose();
            S12.row(i) = imass_psi * (s12_gauss.array() * dcrit.array()).matrix().transpose();
            S22.row(i) = imass_psi * (s22_gauss.array() * dcrit.array()).matrix().transpose();

            // Update damage
            // D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / RefScaleCanada::damage_timescale;

            const Eigen::Matrix<Nextsim::FloatType, 3, 9> imass_psi2 = ParametricTools::massMatrix<3>(smesh, i).inverse()
                * (PSI<3,3>.array().rowwise() * (GAUSSWEIGHTS_3.array() * J.array())).matrix();

            D.row(i) += imass_psi2 * (((1.0 - d_gauss.array()) * (1.0 - dcrit.array()) * dt_mom / RefScaleCanada::damage_timescale).matrix().transpose());

            /*


            //   // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)
            S11.row(i) *= (1.0 - 1.0 / alpha);
            S12.row(i) *= (1.0 - 1.0 / alpha);
            S22.row(i) *= (1.0 - 1.0 / alpha);

            //  zeta = P / 2.0 / DELTA;
            //  eta = zeta / 4;
            //  (zeta-eta) = zeta - 1/4 zeta = 3/4 zeta

            //   S11.row(i)
            //       += 1.0 / alpha * (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S11(i, 0) -= 1.0 / alpha * 0.5 * P;

            const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
            // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // with the gauss weights and with J. This is a 8 x 9 matrix
            const Eigen::Matrix<Nextsim::FloatType, 8, 9> imass_psi = ParametricTools::massMatrix<8>(smesh, i).inverse()
                * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS_3.array() * J.array())).matrix();

            S11.row(i) += imass_psi * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix().transpose());

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += imass_psi * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix().transpose());

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += imass_psi * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix().transpose());


            */
        }
    }

} /* namespace MEB */

} /* namespace Nextsim */

#endif /* __MEB_HPP */
