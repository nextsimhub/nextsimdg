#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "dginitial.hpp"
#include "dgtransport.hpp"
#include "dgvisu.hpp"
#include "mesh.hpp"
#include "newdynamics.hpp"
#include "stopwatch.hpp"

bool WRITE_VTK = true;
int WRITE_EVERY = 10000;

namespace Nextsim {
extern Timer GlobalTimer;
}

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

constexpr double Pstar = 27500.; //!< Ice strength
constexpr double fc = 1.46e-4; //!< Coriolis

constexpr double DeltaMin = 2.e-9; //!< Viscous regime

// parameters form nextsim options.cpp line 302
constexpr double compaction_param = -20; //!< Compation parameter
constexpr double undamaged_time_relaxation_sigma = 1e7; //!< seconds
constexpr double exponent_relaxation_sigma = 5;
constexpr double young = 5.9605e+08;
constexpr double nu0 = 1. / 3.; //!< \param Poisson's ratio
constexpr double compr_strength = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]
constexpr double tan_phi = 0.7; //! \param tan_phi (double) Internal friction coefficient (mu)
constexpr double C_lab = 2.0e6; //! \param C_lab (double) Cohesion at the lab scale (10 cm) [Pa]

}

inline constexpr double SQR(double x) { return x * x; }

//! Description of the problem data, wind & ocean fields
class OceanX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return RefScale::vmax_ocean * (2.0 * y / RefScale::L - 1);
    }
};
class OceanY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return RefScale::vmax_ocean * (1.0 - 2.0 * x / RefScale::L);
    }
};

class AtmX : virtual public Nextsim::InitialBase {
    double time;

public:
    void settime(double t)
    {
        time = t;
    }
    double operator()(double x, double y) const
    {
        constexpr double oneday = 24.0 * 60.0 * 60.0;
        //! Center of cyclone (in m)
        double cM = 256000. + 51200. * time / oneday;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cM) + SQR(y - cM))) * 1.e-3;

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * RefScale::vmax_atm * (cos(alpha) * (x - cM) + sin(alpha) * (y - cM));
    }
};
class AtmY : virtual public Nextsim::InitialBase {
    double time;

public:
    void settime(double t)
    {
        time = t;
    }
    double operator()(double x, double y) const
    {
        constexpr double oneday = 24.0 * 60.0 * 60.0;
        //! Center of cyclone (in m)
        double cM = 256000. + 51200. * time / oneday;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cM) + SQR(y - cM))) * 1.e-3;

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * RefScale::vmax_atm * (-sin(alpha) * (x - cM) + cos(alpha) * (y - cM));
    }
};
class InitialH : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 0.3 + 0.005 * (sin(6.e-5 * x) + sin(3.e-5 * y));
    }
};
class InitialA : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 1.0;
    }
};
class InitialD : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 0.0;
    }
};

int main()
{
    Nextsim::NewDynamics dynamics;

    //! initialize the mesh
    //constexpr size_t N = 25; //!< Number of mesh nodes
    //mesh.BasicInit(N, N, RefScale::L / N, RefScale::L / N);
    //std::cout << "--------------------------------------------" << std::endl;
    //std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    //! Define the spatial mesh
    Nextsim::Mesh mesh;
    constexpr size_t N = 50; //!< Number of mesh nodes
    mesh.BasicInit(N, N, RefScale::L / N, RefScale::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    //! define the time mesh
    const int hours = 48;
    constexpr double T = hours * 60 * 60; //!< Time horizon (in sec) max hours = 2 *  24
    constexpr double dt_adv = 120; // 90 //!< Time step of advection problem
    constexpr size_t NT = T / dt_adv + 1.e-4; //!< Number of Advections steps
    constexpr size_t mom_substeps = 100;
    constexpr double dt_momentum = dt_adv / mom_substeps; //!< Time step of momentum problem

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "Momentym substeps " << mom_substeps << std::endl;

    //! VTK output
    constexpr double T_vtk = 1.0 * 60.0 * 60.0; // evey 1 hours
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 10.0 * 60.0; // every 30 minute
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    // Compute Parabolic CFL
    constexpr double gamma = 1e6; //!< Penalty parameter for internal continuity
    constexpr double gammaboundary = gamma; //!< Penalty parameter for boundary data

    //const double gamma = 5.0; //!< parameter in front of internal penalty terms
    //const double gammaboundary = 5.0; //!< parameter in front of boundary penalty terms

    //! Definition of Variables
    Nextsim::CellVector<2> vx(mesh), vy(mesh); //!< velocity fields

    //! Transport
    Nextsim::DGTransport<2> dgtransport(vx, vy);
    dgtransport.settimesteppingscheme("rk3");
    dgtransport.setmesh(mesh);

    //Damage Transport DGDegree = 0
    //Nextsim::DGTransport<0> dg0transport(vx, vy);
    //dg0transport.settimesteppingscheme("rk1");
    //dg0transport.setmesh(mesh);

    Nextsim::CellVector<1> S11(mesh), S12(mesh), S22(mesh); //!< entries of (symmetric) stress tensor
    Nextsim::CellVector<1> E11(mesh), E12(mesh), E22(mesh); //!< entries of (symmetric) strain stress tensor
    Nextsim::CellVector<2> A(mesh), H(mesh); //!< ice height and ice concentration
    // TODO Damage is DG2 since for transport in must have order of the velocity
    Nextsim::CellVector<2> D(mesh); //!< ice damage. ?? Really dG(0) ??
    //Nextsim::CellVector<2> oceanX(mesh), oceanY(mesh); //!< ocean forcing. ?? Higher order??
    //Nextsim::CellVector<0> atmX(mesh), atmY(mesh); //!< ocean forcing. ?? Higher order??

    //! temporary vectors for time stepping
    Nextsim::CellVector<2> tmpX(mesh), tmpY(mesh);

    //! Initial data of the problem

    Nextsim::L2ProjectInitial(mesh, H, InitialH());
    Nextsim::L2ProjectInitial(mesh, A, InitialA());
    Nextsim::L2ProjectInitial(mesh, D, InitialD());

    Nextsim::CellVector<0> AX(mesh); //!< x-component of atm. velocity
    Nextsim::CellVector<0> AY(mesh); //!< y-component of atm. velocity
    AtmX AtmForcingX; //!< stupid names....
    AtmY AtmForcingY;
    AtmForcingX.settime(0.0);
    AtmForcingY.settime(0.0);
    Nextsim::L2ProjectInitial(mesh, AX, AtmForcingX);
    Nextsim::L2ProjectInitial(mesh, AY, AtmForcingY);

    Nextsim::CellVector<2> OX(mesh); //!< x-component of ocean velocity
    Nextsim::CellVector<2> OY(mesh); //!< y-component of ocean velocity
    Nextsim::L2ProjectInitial(mesh, OX, OceanX());
    Nextsim::L2ProjectInitial(mesh, OY, OceanY());

    //! Initialize the velocity
    vx.zero();
    vy.zero();

    //test only
    //Nextsim::L2ProjectInitial(mesh, vx, AtmForcingX);
    //Nextsim::L2ProjectInitial(mesh, vy, AtmForcingY);

    //save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_dg<2>("Results/vx", 0, vx, mesh);
    Nextsim::VTK::write_dg<2>("Results/vy", 0, vy, mesh);
    Nextsim::VTK::write_dg<1>("Results/S11", 0, S11, mesh);
    Nextsim::VTK::write_dg<1>("Results/S12", 0, S12, mesh);
    Nextsim::VTK::write_dg<1>("Results/S22", 0, S22, mesh);
    Nextsim::VTK::write_dg<2>("Results/A", 0, A, mesh);
    Nextsim::VTK::write_dg<2>("Results/H", 0, H, mesh);
    Nextsim::VTK::write_dg<2>("Results/D", 0, D, mesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    /*
      TIME LOOP
    */
    Nextsim::GlobalTimer.start("time loop");

    for (size_t timestep = 1; timestep <= NT; ++timestep) {
        double time = dt_adv * timestep;
        double timeInMinutes = time / 60.0;
        double timeInHours = time / 60.0 / 60.0;
        double timeInDays = time / 60.0 / 60.0 / 24.;

        if (timestep % NT_log == 0)
            std::cout << "\rAdvection step " << timestep << "\t "
                      << std::setprecision(2)
                      << std::fixed
                      << std::setw(10) << std::right
                      << time << "s\t"
                      << std::setw(8) << std::right
                      << timeInMinutes << "m\t"
                      << std::setw(6) << std::right
                      << timeInHours << "h\t"
                      << std::setw(6) << std::right
                      << timeInDays << "d\t\t" << std::flush;

        Nextsim::GlobalTimer.start("time loop - forcing");
        //! Initial (atm) Forcing (ocean is stationary)
        AtmForcingX.settime(time);
        AtmForcingY.settime(time);
        Nextsim::L2ProjectInitial(mesh, AX, AtmForcingX);
        Nextsim::L2ProjectInitial(mesh, AY, AtmForcingY);

        Nextsim::GlobalTimer.stop("time loop - forcing");

        //! Advection
        Nextsim::GlobalTimer.start("time loop - advection");
        dgtransport.reinitvelocity();
        dgtransport.step(dt_adv, A);
        dgtransport.step(dt_adv, H);
        //dg0transport.step(D);
        dgtransport.step(dt_adv, D);

        A.col(0) = A.col(0).cwiseMin(1.0);
        A.col(0) = A.col(0).cwiseMax(0.0);
        H.col(0) = H.col(0).cwiseMax(0.0);
        Nextsim::GlobalTimer.stop("time loop - advection");

        //! Time step
        Nextsim::GlobalTimer.start("time loop -- mom");
        //! Momentum subcycling
        for (size_t mom_step = 0; mom_step < mom_substeps; ++mom_step) {

            //MOMENTUM EQUATION
            tmpX.zero();
            tmpY.zero();

            // d_t U = ...

            // ocean
            // L/(rho H) * RefScale::F_ocean * |velwater - vel| (velwater-vel)
            tmpX.col(0) += 1 / RefScale::rho_ice * RefScale::F_ocean * ((OX.col(0) - vx.col(0)).array().abs() / H.col(0).array() * OX.col(0).array()).matrix();
            tmpX -= 1 / RefScale::rho_ice * RefScale::F_ocean * (vx.array().colwise() * ((OX.col(0) - vx.col(0)).array().abs() / H.col(0).array())).matrix();
            tmpY.col(0) += 1 / RefScale::rho_ice * RefScale::F_ocean * ((OY.col(0) - vy.col(0)).array().abs() / H.col(0).array() * OY.col(0).array()).matrix();
            tmpY -= 1 / RefScale::rho_ice * RefScale::F_ocean * (vy.array().colwise() * ((OY.col(0) - vy.col(0)).array().abs() / H.col(0).array())).matrix();

            // atm.
            // L/(rho H) * RefScale::F_atm * |velatm| velatm
            tmpX.col(0) += 1 / RefScale::rho_ice * RefScale::F_atm * (AX.col(0).array().abs() / H.col(0).array() * AX.col(0).array()).matrix();
            tmpY.col(0) += 1 / RefScale::rho_ice * RefScale::F_atm * (AY.col(0).array().abs() / H.col(0).array() * AY.col(0).array()).matrix();

            // Coriolis  + surface
            //tmpX.col(0) += RefScale::rho_ice * RefScale::fc * (((OX.col(0) - vx.col(0)).array() * H.col(0).array())).matrix();
            //tmpY.col(0) += RefScale::rho_ice * RefScale::fc * (((OY.col(0) - vy.col(0)).array() * H.col(0).array())).matrix();
            //tmpX += RefScale::rho_ice * RefScale::fc * (((OX - vx).array() * H.array())).matrix();
            //tmpY += RefScale::rho_ice * RefScale::fc * (((OY - vy).array() * H.array())).matrix();

            /* 
             Compute the stress tensor
            */

            // This is already defined in benchmark.cpp
            const double scaleSigma = -1.; //!< -1 since -div(S)  term goes to rhs

            dynamics.addStressTensor(scaleSigma, mesh, tmpX, tmpY, S11, S12, S22); //!<  +div(S, nabla phi) - < Sn, phi>
            dynamics.velocityContinuity(gamma, mesh, tmpX, tmpY, vx, vy); //!< penalize velocity jump in inner edges
            dynamics.velocityDirichletBoundary(gammaboundary, mesh, tmpX, tmpY, vx, vy); //!< no-slip for velocity

            vx += dt_momentum * tmpX;
            vy += dt_momentum * tmpY;

            //DAMAGE EQUATION

            // Sigma = D = sym(grad v)
            dynamics.computeStrainRateTensor(mesh, vx, vy, E11, E12, E22); //!< S = 1/2 (nabla v + nabla v^T)
            //std::cout << "ELOOOO" << std::endl;
            Nextsim::GlobalTimer.start("dyn -- mom -- damage");
#pragma omp parallel for
            for (size_t i = 0; i < mesh.n; ++i) {
                // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
                double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
                //std::cout << Pmax << " " << sigma_n << std::endl;
                double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
                double const time_viscous = RefScale::undamaged_time_relaxation_sigma * std::pow((1. - D(i, 0)) * expC, RefScale::exponent_relaxation_sigma - 1.);

                // Plastic failure tildeP
                double tildeP;
                if (sigma_n < 0.) {
                    //below line copied from nextsim finiteelement.cpp
                    //double const Pmax = std::pow(M_thick[cpt], exponent_compression_factor)*compression_factor*expC;
                    double const Pmax = RefScale::Pstar * pow(H(i, 0), 1.5) * exp(-20.0 * (1.0 - A(i, 0)));
                    // tildeP must be capped at 1 to get an elastic response
                    tildeP = std::min(1., -Pmax / sigma_n);
                } else {
                    tildeP = 0.;
                }

                // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
                // min and - from nextsim
                double const multiplicator = std::min(1. - 1e-12,
                    time_viscous / (time_viscous + dt_momentum * (1. - tildeP)));

                double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

                double const Dunit_factor = 1. / (1. - SQR(RefScale::nu0));
                /* Stiffness matrix
                * 1  nu 0
                * nu 1  0
                * 0  0  (1-nu)
                */
                //M_Dunit[0] = Dunit_factor * 1.;
                //M_Dunit[1] = Dunit_factor * nu0;
                //M_Dunit[3] = Dunit_factor * nu0;
                //M_Dunit[4] = Dunit_factor * 1.;
                //M_Dunit[8] = Dunit_factor * (1. - nu0) ;

                //compute SigmaE
                double SigmaE11 = Dunit_factor * (E11(i, 0) + RefScale::nu0 * E22(i, 0));
                double SigmaE22 = Dunit_factor * (RefScale::nu0 * E11(i, 0) + E22(i, 0));
                double SigmaE12 = 1. / (1 - RefScale::nu0) * E12(i, 0);

                //Elasit prediction Eqn. (32)
                S11(i, 0) += dt_momentum * elasticity * SigmaE11;
                S11(i, 0) *= multiplicator;
                S12(i, 0) += dt_momentum * elasticity * SigmaE12;
                S12(i, 0) *= multiplicator;
                S22(i, 0) += dt_momentum * elasticity * SigmaE22;
                S22(i, 0) *= multiplicator;

                //continiue if stress in inside the failure envelope
                /* Compute the shear and normal stresses, which are two invariants of the internal stress tensor */
                double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
                //update sigma_n
                sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
                //cohesion Eqn. (21)
                //Reference length scale is fixed 0.1 since its cohesion parameter at the lab scale (10 cm)
                double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);
                ; // C_lab;...  : cohesion (Pa)

                // d critical Eqn. (29)
                double dcrit;
                if (sigma_n < -RefScale::compr_strength)
                    dcrit = -RefScale::compr_strength / sigma_n;
                else
                    // M_Cohesion[cpt] depends on local random contribution
                    // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]);
                    // finiteelement.cpp#L3834
                    dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

                /* Calculate the adjusted level of damage */
                if ((0. < dcrit) && (dcrit < 1.)) // sigma_s - tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
                {
                    /* Calculate the characteristic time for damage and damage increment */
                    // M_delta_x[cpt] = mesh.hx ???
                    double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                        / std::sqrt(elasticity);

                    // Eqn. (34)
                    D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                    // Recalculate the new state of stress by relaxing elstically Eqn. (36)
                    S11(i, 0) -= S11(i, 0) * (1. - dcrit) * dt_momentum / td;
                    S12(i, 0) -= S12(i, 0) * (1. - dcrit) * dt_momentum / td;
                    S22(i, 0) -= S22(i, 0) * (1. - dcrit) * dt_momentum / td;
                }
            }
            Nextsim::GlobalTimer.stop("dyn -- mom -- damage");
        }
        Nextsim::GlobalTimer.stop("time loop -- mom");

        //! Output
        if (WRITE_VTK)
            if (timestep % NT_vtk == 0) {
                size_t printstep = timestep / NT_vtk;
                Nextsim::GlobalTimer.start("time loop - i/o");

                Nextsim::VTK::write_dg<2>("Results/vx", printstep, vx, mesh);
                Nextsim::VTK::write_dg<2>("Results/vy", printstep, vy, mesh);
                Nextsim::VTK::write_dg<1>("Results/S11", printstep, S11, mesh);
                Nextsim::VTK::write_dg<1>("Results/S12", printstep, S12, mesh);
                Nextsim::VTK::write_dg<1>("Results/S22", printstep, S22, mesh);
                Nextsim::VTK::write_dg<2>("Results/A", printstep, A, mesh);
                Nextsim::VTK::write_dg<2>("Results/H", printstep, H, mesh);
                Nextsim::VTK::write_dg<2>("Results/D", printstep, D, mesh);
                // Nextsim::VTK::write_dg<0>("Results/ox",printstep,OceanX, mesh);
                // Nextsim::VTK::write_dg<0>("Results/oy",printstep,OceanY, mesh);
                // Nextsim::VTK::write_dg<0>("Results/ax",printstep,AtmX, mesh);
                // Nextsim::VTK::write_dg<0>("Results/ay",printstep,AtmY, mesh);
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
