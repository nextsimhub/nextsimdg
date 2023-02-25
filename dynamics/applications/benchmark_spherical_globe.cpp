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
#include "VectorManipulations.hpp"
#include <map>

bool WRITE_VTK = true;


namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return (x * x); }

namespace ReferenceScale {
  constexpr double R1 = Nextsim::EarthRadius * 0.8; // For z=R1, Dirichlet zero

  constexpr double T = 2.0 * 24 * 60. * 60.; //!< Time horizon 2 days
  constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
  double vmax_atm = 20.0;
  
}



class OceanX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return ReferenceScale::vmax_ocean * (1.0 + 0.5 * sin(2.0*x) * cos(y));
    }
};
class OceanY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean  * cos(3.*y);
    }
};

struct AtmX : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
      return ReferenceScale::vmax_atm * sin(4.0*(x-time/ReferenceScale::T)) *cos(4.0*y);
    }
};
struct AtmY : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
      return -ReferenceScale::vmax_atm * cos(8.0*x)*sin(4.0*y);
    }
};

class InitialH : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return 0.3 + 0.01*sin(4*x)*sin(2*y);
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


//! Creates a rectangular mesh of (nearly) whole earth. Ice for  poles * Ny < iy < (1-poles) * Ny
void create_mesh(Nextsim::ParametricMesh& smesh, size_t Nx, size_t Ny) 
{
  smesh.statuslog = -1;
  smesh.CoordinateSystem = Nextsim::SPHERICAL;
  
  smesh.nx = Nx;
  smesh.ny = Ny;
  smesh.nelements = Nx*Ny;
  smesh.nnodes    = (Nx+1)*(Ny+1);
  smesh.vertices.resize(smesh.nnodes,2);
  
  // z coordinate between -0.8*R and 0.8*R
  const double thetamax = asin(0.8);
  
  for (size_t iy = 0; iy <= Ny; ++iy) // lat 
    for (size_t ix = 0; ix <= Nx; ++ix) // lon
      {
	smesh.vertices(iy*(Nx+1)+ix,0) = -M_PI+2.0 * M_PI*ix/Nx;
	smesh.vertices(iy*(Nx+1)+ix,1) = -thetamax+2.0*thetamax*iy/Ny;
      }

  // ice everywhere
  smesh.landmask.resize(Nx*Ny);
  for (size_t i=0;i<Nx*Ny;++i)
    smesh.landmask[i]=1;

  // dirichlet boundary
  for (auto &it :  smesh.dirichlet)
    it.clear();
  for (size_t i=0;i<Nx;++i)
    {
      smesh.dirichlet[0].push_back(i);
      smesh.dirichlet[2].push_back((Ny-1)*Nx + i);
    }

  // periodic boundary
  smesh.periodic.clear();
  smesh.periodic.resize(1); // 1 segments
  for (size_t i=0;i<Ny;++i)
    smesh.periodic[0].push_back(std::array<size_t,4> ({1, (i+1)*Nx-1, i*Nx, i*(Nx+1)}));
}

template <int CG, int DGadvection>
void run_benchmark(const size_t NX, double distort)
{
  Nextsim::ParametricMesh smesh(Nextsim::SPHERICAL);
    //! Define the spatial mesh
  create_mesh(smesh, 2 * NX, NX);

    //! Compose name of output directory and create it
    std::string resultsdir = "BenchmarkSpherical_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(NX);
    std::filesystem::create_directory(resultsdir);

    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh);

    //! define the time mesh
    constexpr double dt_adv = 120.; //!< Time step of advection problem
    
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 1500.0;
    constexpr double beta = 1500.0;
    constexpr size_t NT_evp = 100; //100;

    //! Rheology-Parameters
    Nextsim::VPParameters VP;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta << std::endl
              << "CG/DG " << CG << "\t" << DGadvection  << std::endl;

    //! VTK output
    constexpr double T_vtk = ReferenceScale::T/48.; // 48.0 * 60.0 * 60.0; // every our
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 10.0 * 60.0; // every 30 minute
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    ////////////////////////////////////////////////// Forcing
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceanx(), OceanX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceany(), OceanY());
    AtmX AtmForcingX;
    AtmY AtmForcingY;
    AtmForcingX.settime(0.0);
    AtmForcingY.settime(0.0);
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AtmForcingX);
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AtmForcingY);


    ////////////////////////////////////////////////// Variables and Initial Values
    Nextsim::DGVector<DGadvection> H(smesh), A(smesh); //!< ice height and concentration
    Nextsim::Interpolations::Function2DG(smesh, H, InitialH());
    Nextsim::Interpolations::Function2DG(smesh, A, InitialA());

    ////////////////////////////////////////////////// i/o of initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    if (1) // write initial?
        if (WRITE_VTK) {
            Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh);
            Nextsim::VTK::write_dg(resultsdir + "/A", 0, A, smesh);
            Nextsim::VTK::write_dg(resultsdir + "/H", 0, H, smesh);
            Nextsim::VTK::write_dg(resultsdir + "/Shear", 0, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
        }
    Nextsim::GlobalTimer.stop("time loop - i/o");

    ////////////////////////////////////////////////// Initialize transport
    Nextsim::DGTransport<DGadvection> dgtransport(smesh);
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

        //////////////////////////////////////////////////
        //! Initialize time-dependent forcing
        Nextsim::GlobalTimer.start("time loop - forcing");
        AtmForcingX.settime(time);
        AtmForcingY.settime(time);
        Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AtmForcingX);
        Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AtmForcingY);
        Nextsim::GlobalTimer.stop("time loop - forcing");

        //////////////////////////////////////////////////
        //! Advection step
        Nextsim::GlobalTimer.start("time loop - advection");

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
                Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", printstep, momentum.GetVx(), momentum.GetVy(), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/A", printstep, A, smesh);
                Nextsim::VTK::write_dg(resultsdir + "/H", printstep, H, smesh);
                Nextsim::VTK::write_dg(resultsdir + "/Delta", printstep, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}

int main()
{
  run_benchmark<1, 3>(128, 0.0);
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

