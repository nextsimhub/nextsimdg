#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "dgvisu.hpp"
#include "dynamics.hpp"
#include "stopwatch.hpp"

#include "benchmark_data.hpp"

bool WRITE_VTK = true;
int WRITE_EVERY = 10000;

namespace Nextsim {
extern Timer GlobalTimer;
}

int main()
{
    Nextsim::Dynamics dynamics;

    //! initialize the mesh
    constexpr size_t N = 100; //!< Number of mesh nodes
    dynamics.GetMesh().BasicInit(N, N, ReferenceScale::L / N, ReferenceScale::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    constexpr double T = 2.0 * 24 * 60 * 60; //!< Time horizon (in sec)
    constexpr double k_adv = 90.0; //!< Time step of advection problem

    // Compute Parabolic CFL
    // get the effective viscosity, e.g. d_t v = 2 eta / rho_ice div(sigma)
    constexpr double eta = 1.e12;

    constexpr double effective_viscosity = 2.0 * eta / ReferenceScale::rho_ice;
    constexpr double gamma = 5; //!< Penalty parameter for internal continuity
    constexpr double gammaboundary = gamma; //!< Penalty parameter for boundary data
    constexpr double cfl = 0.05 / gamma / effective_viscosity;
    //! compute time step. Make the advection step a multiple of it.
    constexpr double k = k_adv / static_cast<int>(0.5 + k_adv / (SQR(ReferenceScale::L / N) * cfl));
    constexpr size_t NT = T / k + 1.e-4;
    constexpr size_t NTadv = k_adv / k + 1.e-4;

    std::cout << "Time step size " << k << "\t" << NT << " time steps" << std::endl
              << "Advection step " << k_adv << "\t every " << NTadv << " momentum steps"
              << std::endl;

    dynamics.GetTimeMesh().BasicInit(NT, k_adv, k);

    //! VTK output every hour
    constexpr double Tvtk = 1.0 * 60.0 * 60.0;
    constexpr size_t NTvtk = Tvtk / k;

    //! Initialize the Dynamical Core (vector sizes, etc.)
    dynamics.BasicInit();

    //! Initial data of the problem
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetH(), InitialH());
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetA(), InitialA());

    //! Initialize the velocity
    dynamics.GetVX().zero();
    dynamics.GetVY().zero();
    //Nextsim::L2ProjectInitial(dynamics.GetMesh(),dynamics.GetS11(), InitialS11());
    //Nextsim::L2ProjectInitial(dynamics.GetMesh(),dynamics.GetS12(), InitialS12());
    //Nextsim::L2ProjectInitial(dynamics.GetMesh(),dynamics.GetS22(), InitialS22());
    dynamics.GetD().zero();

    //! Forcing. Ocean forcing is constant in time.
    AtmX AtmForcingX; //!< stupid names....
    AtmY AtmForcingY;

    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetOceanX(), OceanX());
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetOceanY(), OceanY());

    //dynamics.GetVX().col(0) = dynamics.GetOceanX().col(0);
    //dynamics.GetVY().col(0) = dynamics.GetOceanY().col(0);
    //dynamics.GetVX().col(1) = dynamics.GetOceanX().col(1);
    //dynamics.GetVY().col(1) = dynamics.GetOceanY().col(1);
    //dynamics.GetVX().col(2) = dynamics.GetOceanX().col(2);
    //dynamics.GetVY().col(2) = dynamics.GetOceanY().col(2);

    //save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_dg<2>("Results/vx", 0, dynamics.GetVX(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<2>("Results/vy", 0, dynamics.GetVY(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<1>("Results/S11", 0, dynamics.GetS11(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<1>("Results/S12", 0, dynamics.GetS12(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<1>("Results/S22", 0, dynamics.GetS22(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<2>("Results/A", 0, dynamics.GetA(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<2>("Results/H", 0, dynamics.GetH(), dynamics.GetMesh());
    Nextsim::GlobalTimer.stop("time loop - i/o");

    Nextsim::GlobalTimer.start("time loop");
    size_t advectionstep = 0;

    for (size_t timestep = 1; timestep <= dynamics.GetTimeMesh().N; ++timestep) {
        Nextsim::GlobalTimer.start("time loop - reinit");

        double time = dynamics.GetTimeMesh().dt * timestep; //!< current time in seconds

        if (timestep % 1000000 == 0)
            std::cout << "--- Time step " << timestep << "\t advection step " << advectionstep
                      << "-> day " << time / (24.0 * 60.0 * 60.0) << std::endl;

        // advection step
        if (timestep % NTadv == 0) {
            ++advectionstep;

            //! Initial (atm) Forcing (ocean is stationary)
            AtmForcingX.settime(time);
            AtmForcingY.settime(time);
            Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmX(), AtmForcingX);
            Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmY(), AtmForcingY);
            Nextsim::GlobalTimer.stop("time loop - reinit");

            Nextsim::GlobalTimer.start("dyn");
            Nextsim::GlobalTimer.start("dyn -- adv");
            dynamics.advectionStep();
            Nextsim::GlobalTimer.stop("dyn -- adv");
        }

        //TEST ONLY
        /*
      // sigma = D = sym(grad v)
      dynamics.S11.col(0) = 1. / dynamics.mesh.hx * dynamics.vx.col(1);
      dynamics.S11.col(1) = 1. / dynamics.mesh.hx * 2.*dynamics.vx.col(3);
      dynamics.S11.col(2) = 1. / dynamics.mesh.hx * dynamics.vx.col(5);
      
      dynamics.S12.col(0) = 0.5*(dynamics.vy.col(1)/dynamics.mesh.hx + dynamics.vx.col(2)/dynamics.mesh.hy) ;
      dynamics.S12.col(1) = 0.5 * (dynamics.vx.col(5)/dynamics.mesh.hy + 2.0 * dynamics.vy.col(3)/dynamics.mesh.hx );
      dynamics.S12.col(2) = 0.5 * (2.0 * dynamics.vx.col(4)/dynamics.mesh.hy + dynamics.vy.col(5)/dynamics.mesh.hx );
      
      dynamics.S22.col(0) = 1. / dynamics.mesh.hy * dynamics.vy.col(2);
      dynamics.S22.col(1) = 1. / dynamics.mesh.hy * dynamics.vy.col(5);
      dynamics.S22.col(2) = 1. / dynamics.mesh.hy * 2.*dynamics.vy.col(4);



      Nextsim::VTK::write_dg<1>("Results/S11",1,dynamics.GetS11(), dynamics.GetMesh());
      Nextsim::VTK::write_dg<1>("Results/S12",1,dynamics.GetS12(), dynamics.GetMesh());
      Nextsim::VTK::write_dg<1>("Results/S22",1,dynamics.GetS22(), dynamics.GetMesh());
      Nextsim::VTK::write_dg<2>("Results/vx",1,dynamics.GetVX(), dynamics.GetMesh());
	    Nextsim::VTK::write_dg<2>("Results/vy",1,dynamics.GetVY(), dynamics.GetMesh());
      abort();//end of test
      */

        //! Time step
        //dynamics.step();

        // momentum step
        Nextsim::GlobalTimer.start("dyn -- mom");
        dynamics.momentumSubsteps();
        Nextsim::GlobalTimer.stop("dyn -- mom");
        Nextsim::GlobalTimer.stop("dyn");

        //! Output
        if (WRITE_VTK)
            if (timestep % NTvtk == 0) {
                size_t printstep = timestep / NTvtk;
                Nextsim::GlobalTimer.start("time loop - i/o");

                Nextsim::VTK::write_dg<2>("Results/vx", printstep, dynamics.GetVX(), dynamics.GetMesh());
                Nextsim::VTK::write_dg<2>("Results/vy", printstep, dynamics.GetVY(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<1>("Results/S11", printstep, dynamics.GetS11(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<1>("Results/S12", printstep, dynamics.GetS12(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<1>("Results/S22", printstep, dynamics.GetS22(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<2>("Results/A", printstep, dynamics.GetA(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<2>("Results/H", printstep, dynamics.GetH(), dynamics.GetMesh());

                // Nextsim::VTK::write_dg<0>("Results/D",printstep,dynamics.GetD(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/ox",printstep,dynamics.GetOceanX(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/oy",printstep,dynamics.GetOceanY(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/ax",printstep,dynamics.GetAtmX(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/ay",printstep,dynamics.GetAtmY(), dynamics.GetMesh());
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
