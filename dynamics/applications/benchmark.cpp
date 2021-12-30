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
    size_t N = 25;
    dynamics.GetMesh().BasicInit(N, N, 1. / N, 1. / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    // CFL for Laplace:

    // dt = 0.1 * h^2
    int NT = 1000000;
    double k = 1.0 / N / N * 0.01;

    // //! init time mesh [0 to 2] days
    // float hours = 24.; //24
    // double TMAX = 2.0 * hours * 60.0 * 60.0 / ReferenceScale::T;
    // double k = 10.0 / ReferenceScale::T; //!< time step 10 seconds
    // k = 1. / ReferenceScale::T; //!< This is necessary for Laplace

    // int NT = (static_cast<int>((TMAX / k + 1) / 100 + 1) * 100); //!<  No. of time steps dividable by 100

    //    k = TMAX / NT;

    //dynamics.GetTimeMesh().BasicInit(TMAX,NT,1);
    // call constructor with dt and dt_momentum with k and k^2
    dynamics.GetTimeMesh().BasicInit(NT, k, k);
    //dynamics.GetTimeMesh().BasicInit(NT,k,1./(N*N*4));

    std::cout << "Time mesh of [0," << dynamics.GetTimeMesh().tmax << "] with "
              << dynamics.GetTimeMesh().N << " steps, k = " << dynamics.GetTimeMesh().dt << std::endl;

    double vmax = 0.1 * dynamics.GetMesh().hx / dynamics.GetTimeMesh().dt;
    std::cout << "CFL: maximum ice velocity " << vmax << " (reference) "
              << vmax * ReferenceScale::L / ReferenceScale::T << " (m/s) " << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << std::endl;

    //! Initialize the Dynamical Core (vector sizes, etc.)
    dynamics.BasicInit();

    //! Initial data of the problem
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetH(), InitialH());
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetA(), InitialA());

    //! Initialize the velocity
    dynamics.GetVX().zero();
    dynamics.GetVY().zero();
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetVX(), InitialVX());
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetVY(), InitialVY());

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
    for (size_t timestep = 1; timestep <= dynamics.GetTimeMesh().N; ++timestep) {
        Nextsim::GlobalTimer.start("time loop - reinit");
        double time = dynamics.GetTimeMesh().dt * timestep;
        std::cout << std::endl
                  << "--- Time step " << timestep << "\t"
                  << "-> hour " << time * ReferenceScale::T / (60.0 * 60.0) << std::endl;

        //! Initial (atm) Forcing (ocean is stationary)
        AtmForcingX.settime(time);
        AtmForcingY.settime(time);
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmX(), AtmForcingX);
        Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmY(), AtmForcingY);
        Nextsim::GlobalTimer.stop("time loop - reinit");

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
        dynamics.step();

        //! Output
        if (WRITE_VTK)
            if (timestep % WRITE_EVERY == 0) {
                size_t printstep = timestep / WRITE_EVERY;
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
