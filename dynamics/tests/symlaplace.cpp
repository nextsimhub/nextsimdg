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
class ExactVX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return sin(M_PI * x) * sin(M_PI * y);
    }
};
class ExactVY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y);
    }
};
//! Defines the initial solution
class InitialVX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return sin(M_PI * x) * sin(M_PI * y); // + sin(M_PI * 10.0 * x) * sin(M_PI * 5 * y) * 0.3;
    }
};
class InitialVY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y); // + sin(M_PI * 3. * x) * sin(M_PI * 7 * y) * 0.4;
    }
};

//! Defines the right hand side f = -div(nabla v + nabla v^T)
class FX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        double sx = sin(M_PI * x);
        double sy = sin(M_PI * y);
        double c2x = cos(2. * M_PI * x);
        double c2y = cos(2. * M_PI * y);
        return 3.0 * M_PI * M_PI / 2. * sx * sy - 2.0 * M_PI * M_PI * c2x * c2y;
    }
};
class FY : virtual public Nextsim::InitialBase {

public:
    double operator()(double x, double y) const
    {

        double s2x = sin(2. * M_PI * x);
        double s2y = sin(2. * M_PI * y);
        double cx = cos(M_PI * x);
        double cy = cos(M_PI * y);
        return -M_PI * M_PI / 2. * cx * cy + 6.0 * M_PI * M_PI * s2x * s2y;
    }
};

int main()
{
    Nextsim::Dynamics dynamics;

    //! initialize the mesh
    size_t N = 5;
    double T = 1.0;

    for (int refine = 1; refine <= 3; ++refine) {
        N *= 2;

        double k = 1.0 / N / N * 0.02 / 50.0; // time step size
        size_t NT = static_cast<size_t>(T / k + 1.e-6);

        dynamics.GetMesh().BasicInit(N, N / 1, 1. / N, 1. / N);
        dynamics.GetTimeMesh().BasicInit(NT, k, k);

        WRITE_EVERY = 0.02 / dynamics.GetTimeMesh().dt_momentum;

        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;
        std::cout << "Time mesh with " << NT << " steps, step-size " << k << std::endl;

        //! Initialize the Dynamical Core (vector sizes, etc.)
        dynamics.BasicInit();

        //! Initial data of the problem. First solution
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetVX(), InitialVX());
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetVY(), InitialVY());

        //! right hand side
        Nextsim::CellVector<2> fx(dynamics.GetMesh()), fy(dynamics.GetMesh());
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), fx, FX());
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), fy, FY());

        Nextsim::GlobalTimer.start("time loop");
        //        dynamics.GetTimeMesh().N = 1;

        const double gamma = 50.0; //!< parameter in front of internal penalty terms
        const double gammaboundary = 50; //!< parameter in front of boundary penalty terms

        for (size_t timestep = 1; timestep <= dynamics.GetTimeMesh().N; ++timestep) {

            // v += k * ( F + div( 1/2 (nabla v + nabla v^T) ) )
            dynamics.GetTMPX() = fx;
            dynamics.GetTMPY() = fy;

            dynamics.computeStrainRateTensor(); //!< E = 1/2 (nabla v + nabla v^T)
            dynamics.GetS11() = dynamics.GetE11();
            dynamics.GetS12() = dynamics.GetE12();
            dynamics.GetS22() = dynamics.GetE22();
            dynamics.addStressTensor(-1.0); //!< tmp += div(S)
            dynamics.velocityContinuity(gamma); //!< tmp += < [v], [phi] >
            dynamics.velocityDirichletBoundary(gammaboundary); //!< tmp += < v, [phi] >_G

            dynamics.GetVX() += dynamics.GetTimeMesh().dt_momentum * dynamics.GetTMPX();
            dynamics.GetVY() += dynamics.GetTimeMesh().dt_momentum * dynamics.GetTMPY();

            if (timestep % WRITE_EVERY == 0) {

                double error = sqrt(pow(L2Error(dynamics.GetMesh(), dynamics.GetVX(), ExactVX()), 2.0) + pow(L2Error(dynamics.GetMesh(), dynamics.GetVY(), ExactVY()), 2.0));
                std::cout << "Time step " << timestep << " / " << dynamics.GetTimeMesh().dt_momentum * timestep << std::flush;
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
