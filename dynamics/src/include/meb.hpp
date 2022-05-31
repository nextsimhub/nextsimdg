/*!
 * @file meb.hpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MEB_HPP
#define __MEB_HPP

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

//TODO missing 45\deg it goes to Compresssion
constexpr double tan_phi = 0.7; //! \param tan_phi (double) Internal friction coefficient (mu)
constexpr double sin_phi = 0.7; //! \param sin_phi (double) Internal friction coefficient (mu)

constexpr double c0 = 10.e3; //! \param tan_phi (double) Internal friction coefficient (mu)
constexpr double sigma_c0 = 50.e3; //! \param tan_phi (double) Internal friction coefficient (mu)

//constexpr double C_lab = 2.0e6; //! \param C_lab (double) Test [Pa]
//constexpr double time_relaxation_damage = 2160000.; //!< 25 days in seconds
//
//constexpr double compression_factor = 10e3; //! \param Max pressure for damaged converging ice
//constexpr double exponent_compression_factor
//    = 1.5; //! \param Power of ice thickness in the pressure coefficient

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

    /*
    template <class... CellVectorDGstress, class... CellVectorDGtracer>
    void StressUpdateParameters(const Mesh& mesh, const double dt_momentum,
        CellVectorDGstress... argsStress,
        CellVectorDGtracer... argsTracer)
    {
        //std::cout << argsStress... << std::endl;
        //std::cout << argsTracer... << std::endl;
    }*/

    /**
     * @brief Calculate Stresses for the current time step and update damage.
     * 
     * @details MEB model with additional outputs for development
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
     * @param DELTA delta
     * @param SHEAR shear
     * @param S1 Stress invariant 1
     * @param S2 Stress invariant 2
     * @param eta1 viscosity 1
     * @param eta2 viscosity 2
     * @param stressrelax indicator where stress is relaxed
     * @param sigma_outside indicator where sigma is outside plastic envelope
     * @param tP tilda P
     * @param Pm maximum P
     * @param dt_momentum timestep for momentum subcycle 
     */
    template <int DGstress, int DGtracer>
    void StressUpdateSandbox(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, CellVector<DGtracer>& D, CellVector<DGtracer>& DELTA,
        CellVector<DGtracer>& SHEAR, CellVector<DGtracer>& S1, CellVector<DGtracer>& S2,
        CellVector<DGtracer>& eta1, CellVector<DGtracer>& eta2, CellVector<DGtracer>& stressrelax,
        CellVector<DGtracer>& sigma_outside, CellVector<DGtracer>& tP, CellVector<DGtracer>& Pm,
        CellVector<DGtracer>& Td, CellVector<DGtracer>& d_crit, CellVector<DGtracer>& Regime,
        CellVector<DGtracer>& Multip, CellVector<DGtracer>& Lambda,

        const double dt_momentum)
    {
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Check 1: Not mentioned in the paper
            // There's no ice so we set sigma to 0 and carry on
            const double min_c = 0.1;
            // if (M_conc[cpt] <= min_c) {
            if (A(i, 0) <= min_c) {

                // M_damage[cpt] = 0.;
                D(i, 0) = 0.0;
                // for(int i=0;i<3;i++)
                //     M_sigma[i][cpt] = 0.;
                S11.row(i) *= 0.0;
                S12.row(i) *= 0.0;
                S22.row(i) *= 0.0;
                std::cout << "NO ICE" << std::endl;
                continue;
            }

            //======================================================================
            //! - Updates the internal stress
            //======================================================================

            // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
            double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            // shear stress
            double tau = std::hypot(0.5 * (S11(i, 0) - S22(i, 0)), S12(i, 0));
            assert(tau >= 0);

            const double C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);
            // std::cout << "cfix = " << C_fix << std::endl; std::abort();

            // Outside Envelope
            sigma_outside(i) = 0;
            if (C_fix + RefScale::tan_phi * sigma_n <= tau)
                sigma_outside(i) = tau / (C_fix + RefScale::tan_phi * sigma_n);

            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));

            double const time_viscous = RefScale::undamaged_time_relaxation_sigma
                * std::pow((1. - D(i, 0)) * expC, RefScale::exponent_relaxation_sigma - 1.);

            // Plastic failure tildeP
            double tildeP;
            double const Pmax
                = RefScale::compression_factor * pow(H(i, 0), 1.5) * exp(-20.0 * (1.0 - A(i, 0)));

            // Strange Function make plot
            if (sigma_n < 0.) {
                tildeP = std::min(1., -Pmax / sigma_n);
            } else {
                tildeP = 0.;
            }

            //Region
            if (sigma_n > 0)
                Regime(i) = 3;
            else if (-Pmax > sigma_n)
                Regime(i) = 2;
            else
                Regime(i) = 1;

            // tildeP = 0.0;
            tP(i) = tildeP;
            Pm(i) = Pmax;
            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
            double const multiplicator
                = std::min(1. - 1e-12, time_viscous / (time_viscous + dt_momentum * (1. - tildeP)));

            Multip(i) = multiplicator;
            Lambda(i) = time_viscous;

            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;
            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScale::nu0) * E11.row(i)
                    + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += dt_momentum * elasticity * 1. / (1. + RefScale::nu0) * E12.row(i);
            S22.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScale::nu0) * E22.row(i)
                    + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));

            eta1(i) = elasticity * (1. / (1. + RefScale::nu0));
            eta2(i) = elasticity * Dunit_factor * RefScale::nu0;

            // Eqn. (34)
            S11.row(i) *= multiplicator;
            S12.row(i) *= multiplicator;
            S22.row(i) *= multiplicator;

            //======================================================================
            //! - Estimates the level of damage from the updated internal stress and the local
            //! damage criterion
            //======================================================================

            // continiue if stress in inside the failure envelope
            //  Compute the shear and normal stresses, which are two invariants of the internal
            //  stress tensor
            double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
            // update sigma_n
            sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));

            // Check 3: Discuss Cohesion with Einar
            // Cohesion Eqn. (21)
            // Reference length scale is fixed 0.1 since its cohesion parameter at the lab scale (10
            // cm) C_lab;...  : cohesion (Pa)

            // d critical Eqn. (29)

            double dcrit;
            if (sigma_n < -RefScale::compr_strength)
                dcrit = -RefScale::compr_strength / sigma_n;
            else
                // M_Cohesion[cpt] depends on local random contribution
                // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]);
                // finiteelement.cpp#L3834
                dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            // check minus below
            // dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            // Calculate the characteristic time for damage and damage increment
            // M_delta_x[cpt] = mesh.hx ???
            double td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                / std::sqrt(elasticity);

            Td(i) = td;
            td = 20;
            d_crit(i) = dcrit;
            // Calculate the adjusted level of damage
            stressrelax(i) = 0.0;

            // change that if
            if ((0. < dcrit)
                && (dcrit
                    < 1.)) // sigma_s - tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
            // if (C_fix + RefScale::tan_phi * sigma_n <= sigma_s )
            {
                // std::cout << "DAMAGE" << std::endl;
                // SAVE DAMAGE
                stressrelax(i) = 1.0;

                // Eqn. (34)
                D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                // Recalculate the new state of stress by relaxing elstically Eqn. (36)
                S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
            }

            // Relax damage
            // Check 4: Values are not in the Paper, cf. Eqn. (30)
            // D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage *
            // std::exp(-20. * (1. - A(i, 0))));
            D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage);

            // ellipse-output
            S1(i) = 0.5 * (S11(i, 0) + S22(i, 0));
            S2(i) = std::hypot(S12(i, 0), 0.5 * (S11(i, 0) - S22(i, 0)));
        }
    }
    /*!
     * @brief Stress update for VP Model
     *
     * @details Function used to validate frocing gives same results as mEVP
     * 
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
     * @param Pstar ice strength
     * @param DeltaMin Delta min 2e-9
     * @param dt_momentum timestep for momentum subcycle 
     */
    template <int DGstress, int DGtracer>
    void StressUpdateVP(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, const double Pstar, const double DeltaMin,
        const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            double DELTA = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));

            assert(DELTA > 0);

            //! Ice strength
            double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

            double zeta = P / 2.0 / DELTA;
            double eta = zeta / 4;

            // VP
            S11.row(i) = (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            S11(i, 0) -= 0.5 * P;

            S12.row(i) = (2. * eta * E12.row(i));

            S22.row(i) = (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            S22(i, 0) -= 0.5 * P;
        }
    }

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
            //double const multiplicator
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

            //Calcutate d_critical as Veronique:

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
            //D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage * std::exp(RefScale::compaction_param * (1. - A(i, 0))));
        }
    }

    template <int DGstress, int DGtracer>
    void mMEBStressRelaxation(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, CellVector<DGtracer>& D, CellVector<DGtracer>& A, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

            double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
            // update sigma_n
            double const sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);

            // d critical Eqn. (29)
            double dcrit;
            if (sigma_n < -RefScale::compr_strength)
                dcrit = -RefScale::compr_strength / sigma_n;
            else
                dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                / std::sqrt(elasticity);

            // Calculate the adjusted level of damage
            if ((0. < dcrit) && (dcrit < 1.)) {
                // Eqn. (34)
                D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
            }
            // Relax damage
            D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage);
        }
    }

    template <int DGstress, int DGtracer>
    void mMEBDamageUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, CellVector<DGtracer>& D, CellVector<DGtracer>& A, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

            double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
            // update sigma_n
            double const sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);

            // d critical Eqn. (29)
            double dcrit;
            if (sigma_n < -RefScale::compr_strength)
                dcrit = -RefScale::compr_strength / sigma_n;
            else
                dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                / std::sqrt(elasticity);
            // Calculate the adjusted level of damage
            if ((0. < dcrit) && (dcrit < 1.)) {
                // Eqn. (34)
                D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                //S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                //S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                //S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
            }
            // Relax damage
            //D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage);
        }
    }

    template <int DGstress, int DGtracer>
    void mMEBStressUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, CellVector<DGtracer>& D, const double dt_adv,
        const double alpha, CellVector<DGstress>& S11_mmeb, CellVector<DGstress>& S12_mmeb,
        CellVector<DGstress>& S22_mmeb)
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
                tildeP = std::min(1., -Pmax / sigma_n);
            } else {
                tildeP = 0.;
            }

            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
            //double const multiplicator
            //    = std::min(1. - 1e-12, time_viscous / (time_viscous + dt_adv * (1. - tildeP)));

            //double const multiplicator
            //    = std::min(1. - 1e-12, 1. / (1. + alpha + dt_adv / time_viscous * (1. - tildeP)));

            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            double X = (1. - tildeP) / time_viscous;
            //double Z = (alpha + 1) / dt_adv + X; //alpha + 1. / dt_adv + X; //
            //double Y = alpha + 1. + dt_adv * X;

            /*
            S11.row(i) *= alpha / Y;
            S12.row(i) *= alpha / Y;
            S22.row(i) *= alpha / Y;

            S11.row(i) += 1. / Y * S11_mmeb.row(i) + 1. / Y * dt_adv * elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += 1. / Y * S12_mmeb.row(i) + 1. / Y * dt_adv * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
            S22.row(i) += 1. / Y * S22_mmeb.row(i) + 1. / Y * dt_adv * elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            */

            // Alternative
            //S11.row(i) = 1. / Z * (1. / dt_adv * (alpha * S11.row(i) + S11_mmeb.row(i)) + elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i))));
            //S12.row(i) = 1. / Z * (1. / dt_adv * (alpha * S12.row(i) + S12_mmeb.row(i)) + elasticity * 1 / (1 + RefScale::nu0) * E12.row(i));
            //S22.row(i) = 1. / Z * (1. / dt_adv * (alpha * S22.row(i) + S22_mmeb.row(i)) + elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i))));

            // Step by step
            /*
            S11.row(i) *= (1.0 - 1.0 / alpha);
            S12.row(i) *= (1.0 - 1.0 / alpha);
            S22.row(i) *= (1.0 - 1.0 / alpha);

            S11.row(i) += 1. / alpha * elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += 1. / alpha * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
            S22.row(i) += 1. / alpha * elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            */

            //std::cout << 1 / dt_adv << std::endl;std::exit(0);

            // Step by step
            S11.row(i) *= (1. - (X - 1. / dt_adv) / (alpha + 1. / dt_adv));
            S12.row(i) *= (1. - (X - 1. / dt_adv) / (alpha + 1. / dt_adv));
            S22.row(i) *= (1. - (X - 1. / dt_adv) / (alpha + 1. / dt_adv));

            //add eplitic parts
            S11.row(i) += 1. / (alpha + 1. / dt_adv) * (1. / dt_adv * S11_mmeb.row(i));
            S12.row(i) += 1. / (alpha + 1. / dt_adv) * (1. / dt_adv * S12_mmeb.row(i));
            S22.row(i) += 1. / (alpha + 1. / dt_adv) * (1. / dt_adv * S22_mmeb.row(i));
            //add elasticity
            S11.row(i) += 1. / (alpha + 1. / dt_adv) * elasticity * (1. / (1. + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += 1. / (alpha + 1. / dt_adv) * elasticity * 1. / (1. + RefScale::nu0) * E12.row(i);
            S22.row(i) += 1. / (alpha + 1. / dt_adv) * elasticity * (1. / (1. + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
        }
    }

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

            // Eqn. 25
            double time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma
                * std::pow((1. - D(i, 0)), RefScaleCanada::exponent_relaxation_sigma - 1.) / (H(i, 0) * expC);

            // Eqn. 24 Additional multiplic0cation by H
            double elasticity = RefScaleCanada::young * H(i, 0) * (1. - D(i, 0)) * expC;

            // 1. / (1. + dt / lambda) Eqn. 18
            double const multiplicator = 1. / (1. + dt_momentum / time_viscous);

            time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma;
            elasticity = RefScaleCanada::young;

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScaleCanada::nu0) * E11.row(i)
                    + Dunit_factor * RefScaleCanada::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += dt_momentum * elasticity * 1. / (1. + RefScaleCanada::nu0) * E12.row(i);
            S22.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScaleCanada::nu0) * E22.row(i)
                    + Dunit_factor * RefScaleCanada::nu0 * (E11.row(i) + E22.row(i)));

            // Eqn. (34)
            S12.row(i) *= multiplicator;
            S11.row(i) *= multiplicator;
            S22.row(i) *= multiplicator;

            //======================================================================
            //! - Estimates the level of damage from the updated internal stress and the local
            //! damage criterion
            //======================================================================

            // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
            double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            // shear stress
            double tau = std::sqrt(0.25 * (S11(i, 0) - S22(i, 0)) * (S11(i, 0) - S22(i, 0)) + S12(i, 0) * S12(i, 0));
            assert(tau >= 0);

            //
            const double c = RefScaleCanada::c0 * H(i, 0) * expC;
            const double sigma_c = RefScaleCanada::sigma_c0 * H(i, 0) * expC;

            double dcrit(1.0);
            //Regime(i, 0) = 0;
            //// Mohr-Coulomb condition
            //if (tau + RefScaleCanada::sin_phi * sigma_n - c < 0) {
            //    Regime(i, 0) = 1;
            //}
            //if (sigma_n - tau > sigma_c) {
            //    Regime(i, 0) = 2;
            //}

            if (tau + RefScaleCanada::sin_phi * sigma_n - c > 0)
                dcrit = std::min(1., c / (tau + RefScaleCanada::sin_phi * sigma_n));

            //Relax stress
            S11.row(i) *= dcrit;
            S12.row(i) *= dcrit;
            S22.row(i) *= dcrit;

            //Update damage
            D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / RefScaleCanada::damage_timescale;
        }
    }

    template <int DGstress, int DGtracer>
    void ElasticUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, CellVector<DGtracer>& D, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            //! - Updates the internal stress
            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));

            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;
            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // Elasit prediction Eqn. (32)
            S11.row(i) = RefScale::young
                * (1 / (1 + RefScale::nu0) * E11.row(i)
                    + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) = RefScale::young * 1 / (1 + RefScale::nu0) * E12.row(i);
            S22.row(i) = RefScale::young
                * (1 / (1 + RefScale::nu0) * E22.row(i)
                    + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
        }
    }

    template <int DGstress, int DGtracer>
    void ViscoElasticUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, CellVector<DGtracer>& D, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            double const expC = std::exp(RefScaleCanada::compaction_param * (1. - A(i, 0)));

            // Eqn. 25
            double time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma
                * std::pow((1. - D(i, 0)), RefScaleCanada::exponent_relaxation_sigma - 1.) / (H(i, 0) * expC);

            // Eqn. 24 Additional multiplic0cation by H
            double elasticity = RefScaleCanada::young * H(i, 0) * (1. - D(i, 0)) * expC;

            // 1. / (1. + dt / lambda) Eqn. 18
            double const multiplicator = 1. / (1. + dt_momentum / time_viscous);

            //time_viscous = RefScaleCanada::undamaged_time_relaxation_sigma;
            //elasticity = RefScaleCanada::young;

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScaleCanada::nu0) * E11.row(i)
                    + Dunit_factor * RefScaleCanada::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += dt_momentum * elasticity * 1. / (1. + RefScaleCanada::nu0) * E12.row(i);
            S22.row(i) += dt_momentum * elasticity
                * (1. / (1. + RefScaleCanada::nu0) * E22.row(i)
                    + Dunit_factor * RefScaleCanada::nu0 * (E11.row(i) + E22.row(i)));

            // Eqn. (34)
            S12.row(i) *= multiplicator;
            S11.row(i) *= multiplicator;
            S22.row(i) *= multiplicator;
        }
    }

    //! CG1-DG1
    void ViscoElasticUpdate(const Mesh& mesh, CellVector<3>& S11, CellVector<3>& S12,
        CellVector<3>& S22, const CellVector<3>& E11, const CellVector<3>& E12,
        const CellVector<3>& E22, const CellVector<3>& H,
        const CellVector<3>& A, CellVector<3>& D, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            //!  number of gauss points in each direction equalls 3, for 2d 3^2
            const Eigen::Matrix<double, 1, 9> h_gauss = (H.block<1, 3>(i, 0) * BiG33).array().max(0.0).matrix();
            const Eigen::Matrix<double, 1, 9> a_gauss = (A.block<1, 3>(i, 0) * BiG33).array().max(0.0).min(1.0).matrix();
            const Eigen::Matrix<double, 1, 9> d_gauss = (D.block<1, 3>(i, 0) * BiG33).array().max(0.0).min(1.0).matrix();

            const Eigen::Matrix<double, 1, 9> e11_gauss = E11.block<1, 3>(i, 0) * BiG33;
            const Eigen::Matrix<double, 1, 9> e12_gauss = E12.block<1, 3>(i, 0) * BiG33;
            const Eigen::Matrix<double, 1, 9> e22_gauss = E22.block<1, 3>(i, 0) * BiG33;

            auto expC = (-20.0 * (1.0 - a_gauss.array())).exp();
            // Eqn. 24 Additional multiplic0cation by H
            const Eigen::Matrix<double, 1, 9> elasticity = (RefScaleCanada::young * h_gauss.array() * (1. - d_gauss.array()) * expC.array()).matrix();

            auto powalpha = (1. - d_gauss.array()).pow(RefScaleCanada::exponent_relaxation_sigma - 1.);
            // Eqn. 25
            const Eigen::Matrix<double, 1, 9> time_viscous = (RefScaleCanada::undamaged_time_relaxation_sigma * powalpha.array()).matrix();

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // 1. / (1. + dt / lambda) Eqn. 18
            const Eigen::Matrix<double, 1, 3> multiplicator = (1. / (1. + dt_momentum / time_viscous.array())).matrix() * IBC33;

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e11_gauss.array()).matrix() * IBC33
                + dt_momentum * Dunit_factor * RefScaleCanada::nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array())).matrix() * IBC33;

            S12.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e12_gauss.array()).matrix() * IBC33;

            S22.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e22_gauss.array()).matrix() * IBC33
                + dt_momentum * Dunit_factor * RefScaleCanada::nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array())).matrix() * IBC33;

            // Eqn. (34)
            //! Coefficient-wise multiplication
            S12.row(i).array() *= multiplicator.array();
            S11.row(i).array() *= multiplicator.array();
            S22.row(i).array() *= multiplicator.array();

            //S12.row(i).array() *= 1. / (1. + dt_momentum / RefScaleCanada::undamaged_time_relaxation_sigma);
            //S11.row(i).array() *= 1. / (1. + dt_momentum / RefScaleCanada::undamaged_time_relaxation_sigma);
            //S22.row(i).array() *= 1. / (1. + dt_momentum / RefScaleCanada::undamaged_time_relaxation_sigma);
        }
    }

    //! CG1-DG1
    void StressUpdateMEB(const Mesh& mesh, CellVector<3>& S11, CellVector<3>& S12,
        CellVector<3>& S22, const CellVector<3>& E11, const CellVector<3>& E12,
        const CellVector<3>& E22, const CellVector<3>& H,
        const CellVector<3>& A, CellVector<3>& D, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            //!  number of gauss points in each direction equalls 3, for 2d 3^2
            const Eigen::Matrix<double, 1, 9> h_gauss = (H.block<1, 3>(i, 0) * BiG33).array().max(0.0).matrix();
            const Eigen::Matrix<double, 1, 9> a_gauss = (A.block<1, 3>(i, 0) * BiG33).array().max(0.0).min(1.0).matrix();
            const Eigen::Matrix<double, 1, 9> d_gauss = (D.block<1, 3>(i, 0) * BiG33).array().max(0.0).min(1.0).matrix();

            const Eigen::Matrix<double, 1, 9> e11_gauss = E11.block<1, 3>(i, 0) * BiG33;
            const Eigen::Matrix<double, 1, 9> e12_gauss = E12.block<1, 3>(i, 0) * BiG33;
            const Eigen::Matrix<double, 1, 9> e22_gauss = E22.block<1, 3>(i, 0) * BiG33;

            auto expC = (-20.0 * (1.0 - a_gauss.array())).exp();
            // Eqn. 24 Additional multiplic0cation by H
            const Eigen::Matrix<double, 1, 9> elasticity = (RefScaleCanada::young * h_gauss.array() * (1. - d_gauss.array()) * expC.array()).matrix();

            auto powalpha = (1. - d_gauss.array()).pow(RefScaleCanada::exponent_relaxation_sigma - 1.);
            // Eqn. 25
            const Eigen::Matrix<double, 1, 9> time_viscous = (RefScaleCanada::undamaged_time_relaxation_sigma * powalpha.array()).matrix();

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            // 1. / (1. + dt / lambda) Eqn. 18
            const Eigen::Matrix<double, 1, 3> multiplicator = (1. / (1. + dt_momentum / time_viscous.array())).matrix() * IBC33;

            // Elasit prediction Eqn. (32)
            S11.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e11_gauss.array()).matrix() * IBC33
                + dt_momentum * Dunit_factor * RefScaleCanada::nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array())).matrix() * IBC33;

            S12.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e12_gauss.array()).matrix() * IBC33;

            S22.row(i) += dt_momentum * 1. / (1. + RefScaleCanada::nu0) * (elasticity.array() * e22_gauss.array()).matrix() * IBC33
                + dt_momentum * Dunit_factor * RefScaleCanada::nu0 * (elasticity.array() * (e11_gauss.array() + e22_gauss.array())).matrix() * IBC33;

            // Eqn. (34)
            //! Coefficient-wise multiplication
            S12.row(i).array() *= multiplicator.array();
            S11.row(i).array() *= multiplicator.array();
            S22.row(i).array() *= multiplicator.array();

            //======================================================================
            //! - Estimates the level of damage from the updated internal stress and the local
            //! damage criterion
            //======================================================================

            const Eigen::Matrix<double, 1, 9> s11_gauss = S11.block<1, 3>(i, 0) * BiG33;
            const Eigen::Matrix<double, 1, 9> s12_gauss = S12.block<1, 3>(i, 0) * BiG33;
            const Eigen::Matrix<double, 1, 9> s22_gauss = S22.block<1, 3>(i, 0) * BiG33;

            // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
            //double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            const Eigen::Matrix<double, 1, 9> sigma_n = 0.5 * (s11_gauss.array() + s22_gauss.array());

            // shear stress
            //double tau = std::sqrt(0.25 * (S11(i, 0) - S22(i, 0)) * (S11(i, 0) - S22(i, 0)) + S12(i, 0) * S12(i, 0));

            const Eigen::Matrix<double, 1, 9> tau = (0.25 * (s11_gauss.array() - s22_gauss.array()).square() + s12_gauss.array().square()).sqrt();
            //assert(tau >= 0);

            //const double c = RefScaleCanada::c0 * H(i, 0) * std::exp(RefScaleCanada::compaction_param * (1. - A(i, 0)));
            //const double sigma_c = RefScaleCanada::sigma_c0 * H(i, 0) * std::exp(RefScaleCanada::compaction_param * (1. - A(i, 0)));
            const Eigen::Matrix<double, 1, 9> c = RefScaleCanada::c0 * h_gauss.array() * expC;
            const Eigen::Matrix<double, 1, 9> sigma_c = RefScaleCanada::sigma_c0 * h_gauss.array() * expC;

            /*
            //Regime(i, 0) = 0;
            //// Mohr-Coulomb condition
            //if (tau + RefScaleCanada::sin_phi * sigma_n - c < 0) {
            //    Regime(i, 0) = 1;
            //}
            //if (sigma_n - tau > sigma_c) {
            //    Regime(i, 0) = 2;
            //}
            */
            //double dcrit(1.0);
            Eigen::Matrix<double, 1, 9> dcrit;
            dcrit << 1., 1., 1., 1., 1., 1., 1., 1., 1.;
            //if (tau + RefScaleCanada::sin_phi * sigma_n - c > 0)
            //    dcrit = std::min(1., c / (tau + RefScaleCanada::sin_phi * sigma_n));

            dcrit = (tau.array() + RefScaleCanada::sin_phi * sigma_n.array() - c.array() > 0).select(c.array() / (tau.array() + RefScaleCanada::sin_phi * sigma_n.array()), dcrit);
            auto dcrit_int = dcrit * IBC33;
            //Relax stress
            //S11.row(i) *= dcrit;
            //S12.row(i) *= dcrit;
            //S22.row(i) *= dcrit;
            S11.row(i) = (s11_gauss.array() * dcrit.array()).matrix() * IBC33;
            S12.row(i) = (s12_gauss.array() * dcrit.array()).matrix() * IBC33;
            S22.row(i) = (s22_gauss.array() * dcrit.array()).matrix() * IBC33;

            //Update damage
            //D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / RefScaleCanada::damage_timescale;
            D.row(i) += ((1.0 - d_gauss.array()) * (1.0 - dcrit.array()) * dt_momentum / RefScaleCanada::damage_timescale).matrix() * IBC33;
        }
    }

} /* namespace MEB */

} /* namespace Nextsim */

#endif /* __MEB_HPP */
