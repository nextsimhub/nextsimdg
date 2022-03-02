#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "dginitial.hpp"
#include "dgvisu.hpp"
#include "dynamics.hpp"
#include "stopwatch.hpp"

bool WRITE_VTK = true;
int WRITE_EVERY = 10000;

constexpr double L = 512000.0;
constexpr double eta = 1.e8;
constexpr double rho_ice = 900;

namespace Nextsim {
extern Timer GlobalTimer;
}

/**
 * Solves a heat-like equation
 *
 *   d_t v - div (1/2 (nabla v + nabla v^T) ) = f
 *
 * for testing the convergence of the stress-tensor handling.
 *
 * The exact solution in the stationary limit v(x,y) is given
 * The equation is discretized explicitely in time with very
 * small time steps satisfying k << h^2
 */

//! Defines the exact solution
struct ExactVX {
public:
    double operator()(double x, double y) const { return sin(M_PI * x / L) * sin(M_PI * y / L); }
};
struct ExactVY {
public:
    double operator()(double x, double y) const
    {
        return sin(2.0 * M_PI * x / L) * sin(2.0 * M_PI * y / L);
    }
};
//! Defines the initial solution
struct InitialVX {
public:
    double operator()(double x, double y) const { return sin(M_PI * x / L) * sin(M_PI * y / L); }
};
struct InitialVY {
public:
    double operator()(double x, double y) const
    {
        return sin(2.0 * M_PI * x / L) * sin(2.0 * M_PI * y / L);
    }
};

//! Defines the right hand side f = -div(nabla v + nabla v^T)
struct FX {
public:
    double operator()(double x, double y) const
    {
        double sx = sin(M_PI * x / L);
        double sy = sin(M_PI * y / L);
        double c2x = cos(2. * M_PI * x / L);
        double c2y = cos(2. * M_PI * y / L);
        return (3.0 * M_PI * M_PI / 2. * sx * sy - 2.0 * M_PI * M_PI * c2x * c2y) / L / L;
    }
};
struct FY {

public:
    double operator()(double x, double y) const
    {

        double s2x = sin(2. * M_PI * x / L);
        double s2y = sin(2. * M_PI * y / L);
        double cx = cos(M_PI * x / L);
        double cy = cos(M_PI * y / L);
        return (-M_PI * M_PI / 2. * cx * cy + 6.0 * M_PI * M_PI * s2x * s2y) / L / L;
    }
};

int main()
{
    Nextsim::Dynamics dynamics;

    //! initialize the mesh
    size_t N = 5;
    double T = rho_ice * L * L / 2.0 / eta;
    std::cout << "Time scale:\t" << T << std::endl;

    constexpr double stressscale = 2.0 * eta / rho_ice;
    std::cout << "Stress scale:\t" << stressscale << std::endl;

    constexpr double rhsscale = 2.0 * eta / rho_ice;
    std::cout << "Rhs scale:\t" << rhsscale << std::endl;

    constexpr double exactsolutionnorm = L / sqrt(2.0); //!< ( |u|^2 + |v|^2 )^1/2

    const std::string solver = "EVP";

    for (int refine = 1; refine <= 3; ++refine) {
        N *= 2;

        double h = L / N;

        double gamma = 25.0; //!< parameter in front of internal penalty terms
        double gammaboundary = gamma; //!< parameter in front of boundary penalty terms
        double cfl = 0.05 / gamma;

        double dt = cfl * rho_ice * L * L / 2.0 / eta / N / N;
        double dt_momentum = dt;
        std::cout << "Time step:\t" << dt << std::endl
                  << "Mesh size:\t" << h << std::endl;

        dynamics.GetMesh().BasicInit(N, N, h, h);
        size_t NT = T / dt + 1.e-8;
        WRITE_EVERY = T / 5 / dt_momentum + 1.e-6;

        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;
        std::cout << "Time mesh with " << NT << " steps, step-size " << dt << std::endl;

        //! Initialize the Dynamical Core (vector sizes, etc.)
        dynamics.BasicInit();

        //! Initial data of the problem. First solution
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetVX(), InitialVX());
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetVY(), InitialVY());

        dynamics.GetVX().zero();
        dynamics.GetVY().zero();

        //! right hand side
        Nextsim::CellVector<2> fx(dynamics.GetMesh()), fy(dynamics.GetMesh());
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), fx, FX());
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), fy, FY());

        Nextsim::GlobalTimer.start("time loop");
        //        dynamics.GetTimeMesh().N = 1;

        for (size_t timestep = 1; timestep <= NT; ++timestep) {

            // v += k * ( F + div( 1/2 (nabla v + nabla v^T) ) )
            dynamics.GetTMPX() = fx;
            dynamics.GetTMPY() = fy;
            dynamics.GetTMPX() *= rhsscale;
            dynamics.GetTMPY() *= rhsscale;

            dynamics.computeStrainRateTensor(); //!< E = 1/2 (nabla v + nabla v^T)

            if (solver == "VP") {
                dynamics.GetS11() = dynamics.GetE11();
                dynamics.GetS12() = dynamics.GetE12();
                dynamics.GetS22() = dynamics.GetE22();

                dynamics.addStressTensor(-1.0 * stressscale); //!< tmp += div(S)
                dynamics.velocityContinuity(gamma * stressscale); //!< tmp += < [v], [phi] >
                dynamics.velocityDirichletBoundary(
                    gammaboundary * stressscale); //!< tmp += < v, [phi] >_G
            } else if (solver == "EVP") {
                constexpr double E = 100.0;

                dynamics.GetS11() += E * dt_momentum * stressscale * dynamics.GetE11();
                dynamics.GetS11() *= 1.0 / (1.0 + E * dt_momentum);
                dynamics.GetS12() += E * dt_momentum * stressscale * dynamics.GetE12();
                dynamics.GetS12() *= 1.0 / (1.0 + E * dt_momentum);
                dynamics.GetS22() += E * dt_momentum * stressscale * dynamics.GetE22();
                dynamics.GetS22() *= 1.0 / (1.0 + E * dt_momentum);

                dynamics.addStressTensor(-1.0); //!< tmp += div(S)
                dynamics.velocityContinuity(gamma * stressscale); //!< tmp += < [v], [phi] >
                dynamics.velocityDirichletBoundary(
                    gammaboundary * stressscale); //!< tmp += < v, [phi] >_G

            } else {
                std::cerr << "Solver " << solver << " not known. Only VP / EVP" << std::endl;
                abort();
            }

            dynamics.GetVX() += dt_momentum * dynamics.GetTMPX();
            dynamics.GetVY() += dt_momentum * dynamics.GetTMPY();

            if (timestep % WRITE_EVERY == 0) {

                double error
                    = (sqrt(pow(L2Error(dynamics.GetMesh(), dynamics.GetVX(), ExactVX()), 2.0)
                          + pow(L2Error(dynamics.GetMesh(), dynamics.GetVY(), ExactVY()), 2.0)))
                    / exactsolutionnorm;
                std::cout.precision(10);
                std::cout << "Time step " << timestep << " / " << dt_momentum * timestep
                          << std::flush;
                std::cout << "\terror: " << error << std::endl;
            }
        }
        Nextsim::GlobalTimer.stop("time loop");

        Nextsim::GlobalTimer.start("VTK");
        Nextsim::VTK::write_dg<2>("Results/vx", refine, dynamics.GetVX(), dynamics.GetMesh());
        Nextsim::VTK::write_dg<2>("Results/vy", refine, dynamics.GetVY(), dynamics.GetMesh());

        Nextsim::VTK::write_dg<1>("Results/S11", refine, dynamics.GetS11(), dynamics.GetMesh());
        Nextsim::VTK::write_dg<1>("Results/S12", refine, dynamics.GetS12(), dynamics.GetMesh());
        Nextsim::VTK::write_dg<1>("Results/S22", refine, dynamics.GetS22(), dynamics.GetMesh());
        Nextsim::GlobalTimer.stop("VTK");
        Nextsim::GlobalTimer.print();

        Nextsim::GlobalTimer.reset();
    }
}
