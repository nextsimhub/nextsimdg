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
  
}

class VX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      //      return cos(y)*1.0; // 1 m/s at the equator around z
      return sin(x)*sin(y);
    }
};
class VY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return cos(x);
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
  std::ofstream OUT(meshname.c_str());
  OUT << "ParametricMesh 2.0" << std::endl
      << Nx << "\t" << Ny << std::endl;
  for (size_t iy = 0; iy <= Ny; ++iy) // lat 
    for (size_t ix = 0; ix <= Nx; ++ix) // lon
      OUT << -M_PI+2.0 * M_PI*ix/Nx << "\t" << -M_PI/2.0+M_PI*iy/Ny << std::endl;
  
  OUT << "landmask " << Nx * Ny << std::endl; // no ice on poles :-)
    for (size_t iy = 0; iy < Ny; ++iy) // lat 
      for (size_t ix = 0; ix < Nx; ++ix) // lon
	if ( (iy==0) || (iy==Ny-1) )
	  OUT << 0 << std::endl;
	else
	  OUT << 1 << std::endl;

    // Dirichlet boundary along y0 (bottom) / y1 (top)
    OUT << "dirichlet " << 2*Nx << std::endl; // horizontal
    for (size_t i=0;i<Nx;++i)
      OUT << 1*Nx + i << "\t" << 0 << std::endl; // lower
    for (size_t i=0;i<Nx;++i)
      OUT << (Ny-1)*Nx + i << "\t" << 2 << std::endl; // upper
    
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
    std::string resultsdir = "TestGlobe_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(N);
    std::filesystem::create_directory(resultsdir);
    
    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh, Nextsim::SPHERICAL);
    
    //! define the time mesh
    constexpr double dt_adv = 60*60; //!< Time step of advection problem
    
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

    ////////////////////////////////////////////////// Set Velocity
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetVx(), VX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetVy(), VY());

    ////////////////////////////////////////////////// Variables and Initial Values
    Nextsim::DGVector<DGadvection> H(smesh), A(smesh); //!< ice height and concentration
    Nextsim::Interpolations::Function2DG(smesh, H, InitialH(), Nextsim::SPHERICAL);
    Nextsim::Interpolations::Function2DG(smesh, A, InitialA(), Nextsim::SPHERICAL);


    // (1) Write out initial velocity
    Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh,true);
    Nextsim::VTK::write_dg(resultsdir + "/A", 0, A, smesh,true);
    Nextsim::VTK::write_dg(resultsdir + "/H", 0, H, smesh,true);

    // (2) Project Velocity to DG strain rate tensor
    momentum.ProjectCGVelocityToDGStrain();
    Nextsim::VTK::write_dg(resultsdir + "/E11", 0, momentum.GetE11(), smesh,true);
    Nextsim::VTK::write_dg(resultsdir + "/E12", 0, momentum.GetE12(), smesh,true);
    Nextsim::VTK::write_dg(resultsdir + "/E22", 0, momentum.GetE22(), smesh,true);
    
    abort();
    
    ////////////////////////////////////////////////// Initialize transport
    Nextsim::DGTransport<DGadvection> dgtransport(smesh, Nextsim::SPHERICAL);
    dgtransport.settimesteppingscheme("rk2");

    ////////////////////////////////////////////////// Main Loop
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

        //////////////////////////////////////////////////
        //! Advection step

        // interpolates CG velocity to DG and reinits normal velocity
        dgtransport.prepareAdvection(momentum.GetVx(), momentum.GetVy());

        // performs the transport steps
	dgtransport.step(dt_adv, A);
	dgtransport.step(dt_adv, H);
	
        //! Gauss-point limiting
        Nextsim::LimitMax(A, 1.0);
        Nextsim::LimitMin(A, 0.0);
        Nextsim::LimitMin(H, 0.0);

        //////////////////////////////////////////////////
        momentum.prepareIteration(H, A);
        // MEVP subcycling
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {
	  momentum.mEVPStep(VP, NT_evp, alpha, beta, dt_adv, H, A);
	  //	  momentum.VPStep(VP, dt_adv, H, A);
	  // <- MPI
        }

        //////////////////////////////////////////////////
        if (WRITE_VTK) // Output
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                int printstep = timestep / NT_vtk + 1.e-4;
                Nextsim::VTK::write_cg(resultsdir + "/cgA", printstep, momentum.GetcgA(), smesh,true);
		Nextsim::VTK::write_cg(resultsdir + "/cgH", printstep, momentum.GetcgH(), smesh,true);
		Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", printstep, momentum.GetVx(), momentum.GetVy(), smesh,true);
                Nextsim::VTK::write_dg(resultsdir + "/A", printstep, A, smesh,true);
                Nextsim::VTK::write_dg(resultsdir + "/H", printstep, H, smesh,true);

		Nextsim::VTK::write_dg(resultsdir + "/S11", printstep, momentum.GetS11(), smesh,true);
		Nextsim::VTK::write_dg(resultsdir + "/S12", printstep, momentum.GetS12(), smesh,true);
		Nextsim::VTK::write_dg(resultsdir + "/S22", printstep, momentum.GetS22(), smesh,true);
		
                Nextsim::VTK::write_dg(resultsdir + "/Delta", printstep, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh,true);
                Nextsim::VTK::write_dg(resultsdir + "/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh,true);
		           }
    }

    std::cout << std::endl;
}

int main()
{
  run_benchmark<2, 6>(128);
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
