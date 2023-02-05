/*!
 * @file benchmark_cartesian_globe.cpp
 * @date 24 July 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

/*!
 *
 * Ice dynamics test case the globe in spherical coordinates
 * of the globe. 

 * Ice on [0,2 Pi] * [- 3/8 Pi, 3/8 Pi]
 * Periodic at x=0 and x = 2pi, Dirichlet at y = +/- 10000
 * 
 * Wind is perturbation of (10,0)
 * Ocean (1,0) * (1.0 - (x/ (3/8 Pi)^2) and zero on land
 *
 * Full ice cover, A=1 and H=0.5 
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "DGTransport.hpp"

#include "Tools.hpp"
#include "cgParametricMomentum.hpp"
#include "cgVector.hpp"
#include "dgInitial.hpp"
#include "dgLimit.hpp"
#include "dgVisu.hpp"
#include "mevp.hpp"
#include "stopwatch.hpp"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

bool WRITE_VTK = true;

namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return x * x; }

//////////////////////////////////////////////////// Benchmark testcase from [Mehlmann / Richter, ...]
//! Description of the problem data, wind & ocean fields

namespace ReferenceScale {
constexpr double T = 30.0 * 24.0 * 60.0 * 60.0; //!< Time horizon 1 month
constexpr double VMX = 2.0*M_PI * Nextsim::EarthRadius / T;  // once around
}

class OceanX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return 0.0;
    }
};
class OceanY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return 0.0;
    }
};

struct AtmX : public Nextsim::Interpolations::Function {

public:
  double time;
    double operator()(double x, double y) const
    {
      double mx = -M_PI + 2.0 * M_PI * time / ReferenceScale::T;
      double dx = x-mx;
      if (dx>M_PI)
	dx -= 2.0 * M_PI;
      if (dx<-M_PI)
	dx += 2.0 * M_PI;
      double my = 0.0;
      double scale = exp(-10.0 * (pow(dx,2.0) + pow(y-my,2.0)));
      double alpha = 72.0 / 180.0 * M_PI;
      return ReferenceScale::VMX + scale * 150.0 * (cos(alpha) * dx + sin(alpha) * (y-my));
    }
};
 
struct AtmY : public Nextsim::Interpolations::Function {

public:
  double time;
  double operator()(double x, double y) const
  {
    double mx = -M_PI + 2.0 * M_PI * time / ReferenceScale::T;
    double dx = x-mx;
    if (dx>M_PI)
      dx -= 2.0 * M_PI;
    if (dx<-M_PI)
      dx += 2.0 * M_PI;
    double my = 0.0;
    double scale = exp(-10.0 * (pow(dx,2.0) + pow(y-my,2.0)));
    double alpha = 72.0 / 180.0 * M_PI;
    return scale * 150.0 * (-sin(alpha) * dx + cos(alpha) * (y-my));
  }
};
class InitialH : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {      
      return 0.5 + 0.1 * cos(4.0*x) * sin(4.0*y);
    }
};
class InitialA : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
  {
    return 1.0;
  }
};
//////////////////////////////////////////////////


//! Creates a rectangular mesh of the whole earth. Ice for  poles * Ny < iy < (1-poles) * Ny
void create_mesh(const std::string meshname, size_t Nx, size_t Ny) 
{
  double poles = 1.0/8.0;
    std::ofstream OUT(meshname.c_str());
    OUT << "ParametricMesh 2.0" << std::endl
        << Nx << "\t" << Ny << std::endl;
    for (size_t iy = 0; iy <= Ny; ++iy) // lat 
      for (size_t ix = 0; ix <= Nx; ++ix) // lon
	OUT << -M_PI+2.0 * M_PI*ix/Nx << "\t" << -M_PI/2.0+M_PI*iy/Ny << std::endl;

    OUT << "landmask " << Nx * Ny << std::endl; // no ice on poles :-)
    size_t y0 = poles*Ny;       // first element that is ice
    size_t y1 = (1.0-poles)*Ny; // last element that is ice
    for (size_t iy = 0; iy < Ny; ++iy) // lat 
      for (size_t ix = 0; ix < Nx; ++ix) // lon
	if ( (y0 <= iy) && (iy <= y1) )
	  OUT << 1 << std::endl;
	else
	  OUT << 0 << std::endl;

    // Dirichlet boundary along y0 (bottom) / y1 (top)
    OUT << "dirichlet " << 2*Nx << std::endl; // horizontal
    for (size_t i=0;i<Nx;++i)
      OUT << y0*Nx + i << "\t" << 0 << std::endl; // lower
    for (size_t i=0;i<Nx;++i)
      OUT << y1*Nx + i << "\t" << 2 << std::endl; // upper
    
    OUT << "periodic 1" << std::endl; // Periodic Y-term [1] left/right
    OUT << Ny << std::endl;
    for (size_t i=0;i<Ny;++i)
      OUT << (i+1)*Nx-1 << "\t" << i*Nx << "\t1" << std::endl;

    OUT.close();
}

template <int CG, int DGadvection>
void run_benchmark(const size_t N)
{
    //! Define the spatial mesh
    create_mesh("tmp-spherical.smesh", 2*N, N);
    Nextsim::ParametricMesh smesh(Nextsim::SPHERICAL);
    smesh.readmesh("tmp-spherical.smesh");

    //! Compose name of output directory and create it
    std::string resultsdir = "SphericalGlobe_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(N);
    std::filesystem::create_directory(resultsdir);
    
    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh, Nextsim::SPHERICAL);
    
    //! define the time mesh
    constexpr double dt_adv = 120; //!< Time step of advection problem
    
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 1200.0;
    constexpr double beta = 1200.0;
    constexpr size_t NT_evp = 100;

    //! Rheology-Parameters
    Nextsim::VPParameters VP;
    VP.fc = 0.0;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta << std::endl
              << "CG/DG " << CG << "\t" << DGadvection  << std::endl;

    //! VTK output
    constexpr double T_vtk = ReferenceScale::T/50.;
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 60.0 * 60.0; // every hour
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    ////////////////////////////////////////////////// Forcing
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceanx(), OceanX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceany(), OceanY());
    AtmX AirX;  AirX.time = 0.0;
    AtmY AirY;  AirY.time = 0.0;
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AirX);
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AirY);



    ////////////////////////////////////////////////// Variables and Initial Values
    Nextsim::DGVector<DGadvection> H(smesh), A(smesh); //!< ice height and concentration
    Nextsim::Interpolations::Function2DG(smesh, H, InitialH(), Nextsim::SPHERICAL);
    Nextsim::Interpolations::Function2DG(smesh, A, InitialA(), Nextsim::SPHERICAL);


    
    ////////////////////////////////////////////////// i/o of initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    if (1) // write initial?
        if (WRITE_VTK) {
	  Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh,true);
	    Nextsim::VTK::write_cg_velocity(resultsdir + "/atm", 0, momentum.GetAtmx(), momentum.GetAtmy(), smesh,true);
	    Nextsim::VTK::write_cg_velocity(resultsdir + "/ocean", 0, momentum.GetOceanx(), momentum.GetOceany(), smesh,true);
	    
            Nextsim::VTK::write_dg(resultsdir + "/A", 0, A, smesh,true);
            Nextsim::VTK::write_dg(resultsdir + "/H", 0, H, smesh,true);
            Nextsim::VTK::write_dg(resultsdir + "/Shear", 0, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh,true);
        }
    Nextsim::GlobalTimer.stop("time loop - i/o");

    ////////////////////////////////////////////////// Initialize transport
    Nextsim::DGTransport<DGadvection> dgtransport(smesh, Nextsim::SPHERICAL);
    dgtransport.settimesteppingscheme("rk2");

    ////////////////////////////////////////////////// Main Loop
    Nextsim::GlobalTimer.start("time loop");
    for (size_t timestep = 1; timestep <= NT; ++timestep) {
        double time = dt_adv * timestep;
        double timeInMinutes = time / 60.0;
        double timeInHours = time / 60.0 / 60.0;
        double timeInDays = time / 60.0 / 60.0 / 24.;

        if (timestep % NT_log == 0)
            std::cout << "\rAdvection step " << timestep << "\t " << std::setprecision(2)
                      << std::fixed << std::setw(10) << std::right << time << "s\t" << std::setw(8)
                      << std::right << timeInMinutes << "m\t" << std::setw(6) << std::right
                      << timeInHours << "h\t" << std::setw(6) << std::right << timeInDays << "d\t\t"
                      << std::flush;


	AirX.time = time;
        AirY.time = time;
        Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AirX);
        Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AirY);

	
        //////////////////////////////////////////////////
        //! Advection step
        Nextsim::GlobalTimer.start("time loop - advection");

	//! time depend. forcing
	AirX.time = time;
	AirY.time = time;
	Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AirX);
	Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AirY);

	
        // interpolates CG velocity to DG and reinits normal velocity
        dgtransport.prepareAdvection(momentum.GetVx(), momentum.GetVy());

        // performs the transport steps
	dgtransport.step(dt_adv, A);
	dgtransport.step(dt_adv, H);
	
        //! Gauss-point limiting
        Nextsim::LimitMax(A, 1.0);
        Nextsim::LimitMin(A, 0.0);
        Nextsim::LimitMin(H, 0.0);
        Nextsim::GlobalTimer.stop("time loop - advection");

        //////////////////////////////////////////////////
        Nextsim::GlobalTimer.start("time loop - mevp");
        momentum.prepareIteration(H, A);
        // MEVP subcycling
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {
	  momentum.mEVPStep(VP, NT_evp, alpha, beta, dt_adv, H, A);
	  //	  momentum.VPStep(VP, dt_adv, H, A);
	  // <- MPI
        }
	Nextsim::GlobalTimer.stop("time loop - mevp");

        //////////////////////////////////////////////////
        if (WRITE_VTK) // Output
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                int printstep = timestep / NT_vtk + 1.e-4;
                Nextsim::GlobalTimer.start("time loop - i/o");
		Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", printstep, momentum.GetVx(), momentum.GetVy(), smesh,true);
                Nextsim::VTK::write_dg(resultsdir + "/A", printstep, A, smesh,true);
                Nextsim::VTK::write_dg(resultsdir + "/H", printstep, H, smesh,true);
		
	    // Nextsim::VTK::write_dg(resultsdir + "/S11", 0, momentum.GetS11(), smesh,true);
	    // Nextsim::VTK::write_dg(resultsdir + "/S12", 0, momentum.GetS12(), smesh,true);
	    // Nextsim::VTK::write_dg(resultsdir + "/S22", 0, momentum.GetS22(), smesh,true);
	    // Nextsim::VTK::write_dg(resultsdir + "/E11", 0, momentum.GetE11(), smesh,true);
	    // Nextsim::VTK::write_dg(resultsdir + "/E12", 0, momentum.GetE12(), smesh,true);
	    // Nextsim::VTK::write_dg(resultsdir + "/E22", 0, momentum.GetE22(), smesh,true);

                Nextsim::VTK::write_dg(resultsdir + "/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh,true);
                Nextsim::GlobalTimer.stop("time loop - i/o");
		           }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}

int main()
{
  run_benchmark<2, 1>(1024);
  // int NN[5] = {32,64,128,256,512};
  // for (int n=0;n<5;++n)
  //   {
  //     run_benchmark<1, 1, 3>(NN[n], 0.0);
  //     run_benchmark<1, 3, 3>(NN[n], 0.0);
  //     run_benchmark<1, 6, 3>(NN[n], 0.0);
  //     run_benchmark<2, 1, 8>(NN[n], 0.0);
  //     run_benchmark<2, 3, 8>(NN[n], 0.0);
  //     run_benchmark<2, 6, 8>(NN[n], 0.0);
  //   }
}
