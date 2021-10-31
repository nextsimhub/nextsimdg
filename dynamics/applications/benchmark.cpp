#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>

#include "stopwatch.hpp"
#include "dynamics.hpp"
#include "dgvisu.hpp"

#include "benchmark_data.hpp"

bool WRITE_VTK = true;
int  WRITE_EVERY = 50;


namespace Nextsim
{
  extern Timer GlobalTimer;
}




int main()
{
  Nextsim::Dynamics dynamics;

  //! initialize the mesh
  size_t N = 100;
  dynamics.GetMesh().BasicInit(N,N,1./N,1./N);
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

  //! init time mesh [0 to 2 days
  double TMAX = 2.0 * 24.0 * 60.0 * 60.0 / ReferenceScale::T;
  double    k = 10.0 / ReferenceScale::T; //!< time step 10 seconds
  int NT = (static_cast<int>((TMAX / k + 1) / 100 + 1) * 100); //!<  No time steps dividable by 100
  k = TMAX / NT;
  dynamics.GetTimeMesh().BasicInit(TMAX,NT,1);

  std::cout << "Time mesh of [0," << TMAX << "] with " << NT <<  " steps, k = " << k << std::endl;
  double vmax = 0.1 * dynamics.GetMesh().hx / dynamics.GetTimeMesh().dt;
  std::cout << "CFL: maximum ice velocity " << vmax << " (reference) "
	    << vmax * ReferenceScale::L/ReferenceScale::T << " (m/s) " << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << std::endl;


  //! Initialize the Dynamical Core (vector sizes, etc.)
  dynamics.BasicInit();

  //! Initial data of the problem
  Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetH(), InitialH());
  Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetA(), InitialA());
  dynamics.GetVX().zero();
  dynamics.GetVY().zero();
  dynamics.GetS11().zero();
  dynamics.GetS12().zero();
  dynamics.GetS22().zero();
  dynamics.GetD().zero();

  

  //! Forcing. Ocean forcing is constant in time.
  AtmX AtmForcingX; //!< stupid names....
  AtmY AtmForcingY;
  
  Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetOceanX(), OceanX());
  Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetOceanY(), OceanY());
  

  Nextsim::GlobalTimer.start("time loop");
  for (size_t timestep = 1; timestep <= dynamics.GetTimeMesh().N; ++timestep)
    {
      Nextsim::GlobalTimer.start("time loop - reinit");
      double time = dynamics.GetTimeMesh().dt * timestep;
      std::cout << std::endl << "--- Time step " << timestep << "\t"
		<< "-> hour " << time * ReferenceScale::T/(60.0*60.0) << std::endl;

      //! Initial (atm) Forcing (ocean is stationary)
      AtmForcingX.settime(time);
      AtmForcingY.settime(time);
      Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmX(), AtmForcingX);
      Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmY(), AtmForcingY);
      Nextsim::GlobalTimer.stop("time loop - reinit");      

      //! Time step
      dynamics.step();

      //! Output
      if (WRITE_VTK)
	if (timestep % WRITE_EVERY == 0)
	  {
	    size_t printstep = timestep/WRITE_EVERY;
	    Nextsim::GlobalTimer.start("time loop - i/o");      
	    
	    Nextsim::VTK::write_dg<2>("Results/vx",printstep,dynamics.GetVX(), dynamics.GetMesh());
	    Nextsim::VTK::write_dg<2>("Results/vy",printstep,dynamics.GetVY(), dynamics.GetMesh());
	    // Nextsim::VTK::write_dg<1>("Results/S11",printstep,dynamics.GetS11(), dynamics.GetMesh());
	    // Nextsim::VTK::write_dg<1>("Results/S12",printstep,dynamics.GetS12(), dynamics.GetMesh());
	    // Nextsim::VTK::write_dg<1>("Results/S22",printstep,dynamics.GetS22(), dynamics.GetMesh());
	    Nextsim::VTK::write_dg<2>("Results/A",printstep,dynamics.GetA(), dynamics.GetMesh());
	    Nextsim::VTK::write_dg<2>("Results/H",printstep,dynamics.GetH(), dynamics.GetMesh());
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
