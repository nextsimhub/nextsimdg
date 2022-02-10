/*----------------------------   meb.hpp     ---------------------------*/
#ifndef __meb_HPP
#define __meb_HPP
/*----------------------------   meb.hpp     ---------------------------*/

#include "dgvector.hpp"

namespace RefScale {
// Benchmark testcase from [Mehlmann / Richter, ...]
constexpr double T = 2 * 24. * 60. * 60.; //!< Time horizon 2 days
constexpr double L = 512000.0; //!< Size of domain
constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
constexpr double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind

constexpr double rho_ice = 900.0; //!< Sea ice density
constexpr double rho_atm = 1.3; //!< Air density
constexpr double rho_ocean = 1026.0; //!< Ocean density

constexpr double C_atm = 1.2e-3; //!< Air drag coefficient
constexpr double C_ocean = 5.5e-3; //!< Ocean drag coefficient

constexpr double F_atm = C_atm * rho_atm; //!< effective factor for atm-forcing
constexpr double F_ocean = C_ocean * rho_ocean; //!< effective factor for ocean-forcing

constexpr double Pstar = 27500.0; //!< Ice strength

constexpr double fc = 1.46e-4; //!< Coriolis

constexpr double DeltaMin = 2.e-9; //!< Viscous regime

// parameters form nextsim options.cpp line 302
constexpr double compaction_param = -20.; //!< Compation parameter

//constexpr double undamaged_time_relaxation_sigma = 1e7; //!< seconds
constexpr double undamaged_time_relaxation_sigma = 1e5; //!< Test more viscous

constexpr double time_relaxation_damage = 2160000.; //!< 25 days in seconds

constexpr double compression_factor = 10e3; //! \param Max pressure for damaged converging ice
constexpr double exponent_compression_factor = 1.5; //! \param Power of ice thickness in the pressure coefficient


constexpr double exponent_relaxation_sigma = 5;
constexpr double young = 5.9605e+08; 
constexpr double nu0 = 1. / 3.; //!< \param Poisson's ratio
constexpr double compr_strength = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]
constexpr double tan_phi = 0.7; //! \param tan_phi (double) Internal friction coefficient (mu)
//constexpr double C_lab = 2.0e6; //! \param C_lab (double) Cohesion at the lab scale (10 cm) [Pa]
constexpr double C_lab = 2.0e6; //! \param C_lab (double) Test [Pa]

}

namespace Nextsim {

/*!
   * This namespace collects the routines required for the MEB solver
   */
namespace MEB {

    inline constexpr double SQR(double x)
    {
        return x * x;
    }

    template <int DGstress, int DGtracer>
    void StressUpdateSandbox(const Mesh& mesh,
        CellVector<DGstress>& S11, CellVector<DGstress>& S12, CellVector<DGstress>& S22,
        const CellVector<DGstress>& E11, const CellVector<DGstress>& E12, const CellVector<DGstress>& E22,
        const CellVector<DGtracer>& H, const CellVector<DGtracer>& A, CellVector<DGtracer>& D,
        CellVector<DGtracer>& DELTA, CellVector<DGtracer>& SHEAR, CellVector<DGtracer>& S1, 
        CellVector<DGtracer>& S2, CellVector<DGtracer>& eta1, CellVector<DGtracer>& eta2,
        CellVector<DGtracer>& stressrelax, CellVector<DGtracer>& sigma_outside,
        CellVector<DGtracer>& tP, CellVector<DGtracer>& Pm,
        const double dt_momentum)
    {
#pragma omp parallel for
            for (size_t i = 0; i < mesh.n; ++i) {

                // Check 1: Not mentioned in the paper
                // There's no ice so we set sigma to 0 and carry on
                const double min_c = 0.1;
                //if (M_conc[cpt] <= min_c) {
                if (A(i, 0) <= min_c) {

                    //M_damage[cpt] = 0.;
                    D(i, 0) = 0.0;

                    //for(int i=0;i<3;i++)
                    //    M_sigma[i][cpt] = 0.;
                    S11.row(i) *= 0.0;
                    S12.row(i) *= 0.0;
                    S22.row(i) *= 0.0;
                    std::cout << "NO ICE" << std::endl;
                    continue;
                }
                
                /*
                // For comparison with mEVP
                DELTA(i, 0) = sqrt(
                    SQR(RefScale::DeltaMin)
                    + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                    + 1.50 * E11(i, 0) * E22(i, 0)
                    + SQR(E12(i, 0)));
                if (!(DELTA(i, 0) > 0)){
                    std::cout << DELTA(i, 0) << std::endl;
                    assert(DELTA(i, 0) > 0);
                }

                //! Ice strength
                double P = RefScale::Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
                SHEAR(i, 0) = sqrt((SQR(E11(i, 0) - E22(i, 0)) + 4.0 * SQR(E12(i, 0))));
                //double zeta = P / 2.0 / DELTA(i, 0);
                //double eta = zeta / 4;
                // save mu and nu fro mEVP
                //mu1(i, 0) = 2*eta ;  mu2(i, 0) = zeta - eta;
                */

                //======================================================================
                //! - Updates the internal stress
                //======================================================================
                
                // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
                double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
                //shear stress
                double tau = std::sqrt( SQR(0.5 * (S11(i, 0) - S22(i, 0))) + SQR(S12(i,0))  );
                // /assert(tau > 0);
                const double C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);

                //std::cout << "cfix = " << C_fix << std::endl; std::abort();


                // Outside Envelope
                sigma_outside(i) = 0;
                if (C_fix + RefScale::tan_phi * sigma_n <= tau )
                  sigma_outside(i) = tau/(C_fix + RefScale::tan_phi * sigma_n);
                  

                double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));

                // Check 1
                double const time_viscous = RefScale::undamaged_time_relaxation_sigma * std::pow((1. - D(i, 0)) * expC, RefScale::exponent_relaxation_sigma - 1.);
                //double const time_viscous = RefScale::undamaged_time_relaxation_sigma * std::pow(H(i, 0) * expC, RefScale::exponent_relaxation_sigma - 1.);

                // Plastic failure tildeP
                double tildeP;
                double const Pmax = RefScale::compression_factor * pow(H(i, 0), 1.5) * exp(-20.0 * (1.0 - A(i, 0)));
                if (sigma_n < 0.) {
                    
                    //double const Pmax = RefScale::Pstar * pow(H(i, 0), 1.5) * exp(-20.0 * (1.0 - A(i, 0)));
                    //double const Pmax = RefScale::compression_factor * pow(H(i, 0), RefScale::exponent_compression_factor) * exp(RefScale::compaction_param * (1.0 - A(i, 0)));
                    // skipping 1.5 power
                    //double const Pmax = RefScale::Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

                    // tildeP must be capped at 1 to get an elastic response
                    tildeP = std::min(1., -Pmax / sigma_n);
                } else {
                    tildeP = 0.;
                }

                //tildeP = 0.0;

                tP(i) = tildeP;
                Pm(i) = Pmax;
                // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
                // Check 2: min 1.-1e-12 Not mentioned in the paper from nextsim
                double const multiplicator = std::min(1. - 2e-4,
                    time_viscous / (time_viscous + dt_momentum * (1. - tildeP)));
                // -12 -> -3

                
                //if ( time_viscous < 100000-1)
                //std::cout << "ELO1" << time_viscous / (time_viscous + dt_momentum * (1. - tildeP)) << std::endl;
                //std::cout << "" << time_viscous / (time_viscous + dt_momentum * (1. - tildeP)) << std::endl;

                // NEW 1 Initially elastcity is deformed by H
                double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;
                //double const elasticity = RefScale::young * H(i, 0) * expC;

                double const Dunit_factor = 1. / (1. - SQR(RefScale::nu0));


                //Elasit prediction Eqn. (32)
                S11.row(i) += dt_momentum * elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
                S12.row(i) += dt_momentum * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
                S22.row(i) += dt_momentum * elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
                
                eta1(i) = elasticity * (1 / (1 + RefScale::nu0));
                eta2(i) = elasticity * Dunit_factor * RefScale::nu0;

                S11.row(i) *= multiplicator;
                S12.row(i) *= multiplicator;
                S22.row(i) *= multiplicator;

                //======================================================================
                //! - Estimates the level of damage from the updated internal stress and the local damage criterion
                //======================================================================
                 
                //continiue if stress in inside the failure envelope
                // Compute the shear and normal stresses, which are two invariants of the internal stress tensor
                double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
                //update sigma_n
                sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));

                // Check 3: Discuss Cohesion with Einar
                //Cohesion Eqn. (21)
                //Reference length scale is fixed 0.1 since its cohesion parameter at the lab scale (10 cm)
                // C_lab;...  : cohesion (Pa)

                // d critical Eqn. (29)

                double dcrit;
                if (sigma_n < -RefScale::compr_strength)
                    dcrit = -RefScale::compr_strength / sigma_n;
                else
                    // M_Cohesion[cpt] depends on local random contribution
                    // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]);
                    // finiteelement.cpp#L3834
                    dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);
                    
                    //check minus below
                    //dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

                // Calculate the characteristic time for damage and damage increment
                // M_delta_x[cpt] = mesh.hx ???
                double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                    / std::sqrt(elasticity);
                
                // Calculate the adjusted level of damage
                stressrelax(i) = 0.0;
                
                //change that if
                if ((0. < dcrit) && (dcrit < 1.)) // sigma_s - tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
                //if (C_fix + RefScale::tan_phi * sigma_n <= sigma_s )
                {
                    //std::cout << "DAMAGE" << std::endl;
                    //SAVE DAMAGE
                    stressrelax(i) = 1.0;

                    // Eqn. (34)
                    D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                    // Recalculate the new state of stress by relaxing elstically Eqn. (36)
                    S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                    S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                    S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
                }

                //Relax damage
                //Check 4: Values are not in the Paper, cf. Eqn. (30)
                //D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage * std::exp(-20. * (1. - A(i, 0))));
                D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage );
            
                // ellipse-output
                S1(i) = 0.5*(S11(i, 0) + S22(i, 0)) ;
                S2(i) = sqrt(SQR(S12(i, 0)) + SQR(0.5 * (S11(i, 0) - S22(i, 0)))) ;
                
            }

    }


    template <int DGstress, int DGtracer>
    void StressUpdateVP(const Mesh& mesh,
        CellVector<DGstress>& S11, CellVector<DGstress>& S12, CellVector<DGstress>& S22,
        const CellVector<DGstress>& E11, const CellVector<DGstress>& E12, const CellVector<DGstress>& E22,
        const CellVector<DGtracer>& H, const CellVector<DGtracer>& A,
        const double Pstar, const double DeltaMin,
        const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            double DELTA = sqrt(
                SQR(DeltaMin)
                + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0)
                + SQR(E12(i, 0)));
            assert(DELTA > 0);

            //! Ice strength
            double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

            double zeta = P / 2.0 / DELTA;
            double eta = zeta / 4;


            //VP
            S11.row(i) =  (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            S11(i, 0) -=  0.5 * P;

            S12.row(i) = (2. * eta * E12.row(i));

            S22.row(i) =  (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            S22(i, 0) -=  0.5 * P;

        }
    }




    template <int DGstress, int DGtracer>
    void StressUpdate(const Mesh& mesh,
        CellVector<DGstress>& S11, CellVector<DGstress>& S12, CellVector<DGstress>& S22,
        const CellVector<DGstress>& E11, const CellVector<DGstress>& E12, const CellVector<DGstress>& E22,
        const CellVector<DGtracer>& H, const CellVector<DGtracer>& A, CellVector<DGtracer>& D,
        const double dt_momentum)
    {

            //! Stress Update
#pragma omp parallel for
            for (size_t i = 0; i < mesh.n; ++i) {

                // Check 1: Not mentioned in the paper
                // There's no ice so we set sigma to 0 and carry on
                const double min_c = 0.1;
                //if (M_conc[cpt] <= min_c) {
                if (A(i, 0) <= min_c) {

                    //M_damage[cpt] = 0.;
                    D(i, 0) = 0.0;

                    //for(int i=0;i<3;i++)
                    //    M_sigma[i][cpt] = 0.;
                    S11.row(i) *= 0.0;
                    S12.row(i) *= 0.0;
                    S22.row(i) *= 0.0;
                    continue;
                }


                /*======================================================================
                //! - Updates the internal stress
                *======================================================================
                */

                // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
                double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
                double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
                double const time_viscous = RefScale::undamaged_time_relaxation_sigma * std::pow((1. - D(i, 0)) * expC, RefScale::exponent_relaxation_sigma - 1.);
             
                // Plastic failure tildeP
                double tildeP;
                if (sigma_n < 0.) {
                    double const Pmax = RefScale::compression_factor * pow(H(i, 0), RefScale::exponent_compression_factor) * exp(RefScale::compaction_param * (1.0 - A(i, 0)));
             
                    // tildeP must be capped at 1 to get an elastic response
                    tildeP = std::min(1., -Pmax / sigma_n);
                } else {
                    tildeP = 0.;
                }

                // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
                // Check 2: min 1.-1e-12 Not mentioned in the paper from nextsim
                double const multiplicator = std::min(1. - 1e-12,
                    time_viscous / (time_viscous + dt_momentum * (1. - tildeP)));
                double const elasticity = RefScale::young * (1. - D(i, 0) ) * expC;

                double const Dunit_factor = 1. / (1. - SQR(RefScale::nu0));


                //Elasit prediction Eqn. (32)
                S11.row(i) += dt_momentum * elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
                S12.row(i) += dt_momentum * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
                S22.row(i) += dt_momentum * elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
                S11.row(i) *= multiplicator;
                S12.row(i) *= multiplicator;
                S22.row(i) *= multiplicator;

                /*======================================================================
                 //! - Estimates the level of damage from the updated internal stress and the local damage criterion
                 *======================================================================
                 */

                //continiue if stress in inside the failure envelope
                // Compute the shear and normal stresses, which are two invariants of the internal stress tensor
                double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
                //update sigma_n
                sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));

                // Check 3: Discuss Cohesion with Einar
                //Cohesion Eqn. (21)
                //Reference length scale is fixed 0.1 since its cohesion parameter at the lab scale (10 cm)
                double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);
                // C_lab;...  : cohesion (Pa)

                // d critical Eqn. (29)
                double dcrit;
                if (sigma_n < -RefScale::compr_strength)
                    dcrit = -RefScale::compr_strength / sigma_n;
                else
                    // M_Cohesion[cpt] depends on local random contribution
                    // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]);
                    // finiteelement.cpp#L3834
                    dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

                // Calculate the characteristic time for damage and damage increment
                // M_delta_x[cpt] = mesh.hx ???
                double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                    / std::sqrt(elasticity);
                // Calculate the adjusted level of damage
                if ((0. < dcrit) && (dcrit < 1.)) // sigma_s - tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
                {
                    // Eqn. (34)
                    D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                    // Recalculate the new state of stress by relaxing elstically Eqn. (36)
                    S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                    S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                    S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
                }

                //Relax damage
                //Check 4: Values are not in the Paper, cf. Eqn. (30)
                D(i, 0) = std::max(0., D(i, 0) - dt_momentum / 2160000. * std::exp(RefScale::compaction_param * (1. - A(i, 0))) );

            }
}




}
}

/*----------------------------   meb.hpp     ---------------------------*/
/* end of #ifndef __men_HPP */
#endif
/*----------------------------   meb.hpp     ---------------------------*/
