/*!
 * @file benchmark_mehlmann_mevp.cpp
 * @date 24 July 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
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

Nextsim::COORDINATES CoordinateSystem = Nextsim::CARTESIAN;

namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return x * x; }

//////////////////////////////////////////////////// Benchmark testcase from [Mehlmann / Richter, ...]
//! Description of the problem data, wind & ocean fields

namespace ReferenceScale {
constexpr double T = 8.0 * 24 * 60. * 60.; //!< Time horizon 8 days
  
constexpr double L = 512000.0; //!< Size of domain !!!
constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind
}

class OceanX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (2.0 * y / ReferenceScale::L - 1.0);
    }
};
class OceanY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (1.0 - 2.0 * x / ReferenceScale::L);
    }
};

struct AtmX : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
        constexpr double oneday = 24.0 * 60.0 * 60.0;
        //! Center of cyclone (in m)
        double cMx = 200000. + 150000. * sin(2.0*M_PI*time/ReferenceScale::T);
	double cMy = 350000. - 150000. * cos(2.0*M_PI*time/ReferenceScale::T);
        double alpha = (90.0 + 18.0 * tanh(20.0*(time/ReferenceScale::T-0.5))) / 180.0 * M_PI;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cMx) + SQR(y - cMy))) * 1.e-3;


	
        return -scale * ReferenceScale::vmax_atm * (cos(alpha) * (x - cMx) + sin(alpha) * (y - cMy));
    }
};
struct AtmY : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
        constexpr double oneday = 24.0 * 60.0 * 60.0;
        //! Center of cyclone (in m)
        //! Center of cyclone (in m)
        double cMx = 200000. + 150000. * sin(2.0*M_PI*time/ReferenceScale::T);
	double cMy = 350000. - 150000. * cos(2.0*M_PI*time/ReferenceScale::T);
        double alpha = (90.0 + 18.0 * tanh(20.0*(time/ReferenceScale::T-0.5))) / 180.0 * M_PI;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cMx) + SQR(y - cMy))) * 1.e-3;

        return -scale * ReferenceScale::vmax_atm * (-sin(alpha) * (x - cMx) + cos(alpha) * (y - cMy));
    }
};
class InitialH : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return 0.3 + 0.005 * (sin(6.e-5 * x) + sin(3.e-5 * y));
    }
};
class InitialA : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const { return 1.0; }
};

template <int CG, int DGadvection>
void run_benchmark()
{
  Nextsim::ParametricMesh smesh(CoordinateSystem);
  smesh.readmesh("../tests/example3-25km_NH.smesh");
  // transform mesh to 2d-projection 
  double R = 6371000.0; 
  for (size_t i = 0; i<smesh.nnodes;++i)
    {
      const double x = R * cos(M_PI/180.0*smesh.vertices(i,1)) * cos(M_PI/180.0*smesh.vertices(i,0));
      const double y = R * cos(M_PI/180.0*smesh.vertices(i,1)) * sin(M_PI/180.0*smesh.vertices(i,0));
      smesh.vertices(i,0) =0.25*(-sin(M_PI/4.)*x+cos(M_PI/4.0)*y       );
      smesh.vertices(i,1) =0.25*( cos(M_PI/4.)*x+sin(M_PI/4.0)*y + 2.e6);
    }

  // output land mask
  Nextsim::DGVector<1> landmask(smesh);
  for (size_t i=0;i<smesh.nelements;++i)
    landmask(i,0) = smesh.landmask[i];
  Nextsim::VTK::write_dg<1>("landmask", 0, landmask, smesh);
				   
  
  // output boundary info
  Nextsim::DGVector<1> boundary(smesh);
  for (size_t j=0;j<4;++j)
    for (size_t i=0;i<smesh.dirichlet[j].size();++i)
      boundary(smesh.dirichlet[j][i],0) = 1+j;
  Nextsim::VTK::write_dg<1>("boundary", 0, boundary, smesh);
  
  
  //! Compose name of output directory and create it
  std::string resultsdir = "BenchmarkOrca_" + std::to_string(CG) + "_" + std::to_string(DGadvection);
  std::filesystem::create_directory(resultsdir);
  
  //! Main class to handle the momentum equation. This class also stores the CG velocity vector
  Nextsim::CGParametricMomentum<CG> momentum(smesh);
  
  //! define the time mesh
  constexpr double dt_adv = 60; //!< Time step of advection problem
  
  constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps
  
  //! MEVP parameters
  constexpr double alpha = 1500.0;
  constexpr double beta = 1500.0;
  constexpr size_t NT_evp = 100;
  
  //! Rheology-Parameters
  Nextsim::VPParameters VP;
  
  std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
	    << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta << std::endl
	    << "CG/DG " << CG << "\t" << DGadvection  << std::endl;
  
  //! VTK output
  constexpr double T_vtk = 0.01*ReferenceScale::T; // 48.0 * 60.0 * 60.0; // once
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
	    Nextsim::VTK::write_cg_velocity(resultsdir + "/atm", 0, momentum.GetAtmx(), momentum.GetAtmy(), smesh);
	    Nextsim::VTK::write_cg_velocity(resultsdir + "/ocean", 0, momentum.GetOceanx(), momentum.GetOceany(), smesh);
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

	//        if (timestep % NT_log == 0)
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
                Nextsim::VTK::write_dg(resultsdir + "/E11", printstep, momentum.GetE11(), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/E12", printstep, momentum.GetE12(), smesh);
		Nextsim::VTK::write_dg(resultsdir + "/E22", printstep, momentum.GetE22(), smesh);
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
  run_benchmark<1, 1>();
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
