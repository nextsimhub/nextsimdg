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

/*!
 *
 * The benchmark [Mehlmann, 2017, 2021] on the spherical orca NH 25km mesh
 * 
 * Mesh rotated to Greenland (Pi/2 lat is in Greenland), pole is at (lat=75^o, lon=0)
 * Pole-Vector is PV = R (cos(75o) cos(40o), cos(75o)sin(40o), sin(75o) )
 *
 * The benchmark data is first transformed from [0,512km]^2 to [-256,256km]^2
 * and then mapped to the domain as plane orthogonal to PV. 
 * 
 * This plane is spanned by the vectors
 * (0,1,0) and (-sin(75o),0,cos(75o) )
 * 
 */

constexpr double L0 = 75./180.*M_PI;
constexpr double L1 = -180.0/180.*M_PI;

constexpr double P1 = Nextsim::EarthRadius * cos(L0) * cos(L1);
constexpr double P2 = Nextsim::EarthRadius * cos(L0) * sin(L1);
constexpr double P3 = Nextsim::EarthRadius * sin(L0);

constexpr double PA1 = -cos(L0) * sin(L1);
constexpr double PA2 =  cos(L0) * cos(L1);
constexpr double PA3 =  0.0;

constexpr double PB1 = (P2*PA3-P3*PA2) / Nextsim::EarthRadius;
constexpr double PB2 = (P3*PA1-P1*PA3) / Nextsim::EarthRadius;
constexpr double PB3 = (P1*PA2-P2*PA1) / Nextsim::EarthRadius;


Nextsim::COORDINATES CoordinateSystem = Nextsim::SPHERICAL;

namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return x * x; }

//////////////////////////////////////////////////// Benchmark testcase from [Mehlmann / Richter, ...]
//! Description of the problem data, wind & ocean fields

namespace ReferenceScale {
constexpr double T = 2.0 * 24 * 60. * 60.; //!< Time horizon 8 days
  
constexpr double L = 512000.0; //!< Size of domain !!!
constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind
}

class OceanX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      double X = Nextsim::EarthRadius * cos(x)*cos(y);
      double Y = Nextsim::EarthRadius * sin(x)*cos(y);
      double Z = Nextsim::EarthRadius * sin(y);
      double xA = X * PA1 + Y * PA2 + Z * PA3;
      double xB = X * PB1 + Y * PB2 + Z * PB3;

      double oA = ReferenceScale::vmax_ocean * 2.0 * xB / ReferenceScale::L;
      double oB = ReferenceScale::vmax_ocean * (-2.0) * xA / ReferenceScale::L;
      // vel points in ox (0,-1,0) + oy(-sin,0,cos) - direction. and must be mapped to lat/lon system

      double ex1 = -sin(x); // elon
      double ex2 =  cos(x);

      return oA * (ex1*PA1 + ex2 * PA2) + oB * (ex1 * PB1 + ex2 * PB2);
    }
};
class OceanY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      double X = Nextsim::EarthRadius * cos(x)*cos(y);
      double Y = Nextsim::EarthRadius * sin(x)*cos(y);
      double Z = Nextsim::EarthRadius * sin(y);
      double xA = X * PA1 + Y * PA2 + Z * PA3;
      double xB = X * PB1 + Y * PB2 + Z * PB3;

      double oA = ReferenceScale::vmax_ocean * 2.0 * xB / ReferenceScale::L;
      double oB = ReferenceScale::vmax_ocean * (-2.0) * xA / ReferenceScale::L;
      // vel points in ox (0,-1,0) + oy(-sin,0,cos) - direction. and must be mapped to lat/lon system

      double ey1 = -sin(y)*cos(x);
      double ey2 = -sin(y)*sin(x);
      double ey3 = cos(y);

      return oA * (ey1*PA1 + ey2 * PA2 + ey3 * PA3) + oB * (ey1 * PB1 + ey2 * PB2 + ey3 * PB3);
      
    }
};

struct AtmX : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
      // xA,xB coordinate in Arctic plane
      double X = Nextsim::EarthRadius * cos(x)*cos(y);
      double Y = Nextsim::EarthRadius * sin(x)*cos(y);
      double Z = Nextsim::EarthRadius * sin(y);
      double xA = X * PA1 + Y * PA2 + Z * PA3;
      double xB = X * PB1 + Y * PB2 + Z * PB3;
      

      constexpr double oneday = 24.0 * 60.0 * 60.0;
      //! Center of cyclone (in m) (cmX, cmX
      double cMx = 51200.0 * time/oneday;
      double cMxA = cMx * PA1 + cMx * PA2; // same in the plane-system
      double cMxB = cMx * PA2 + cMx * PA2;

      
      
      double alpha = 72.0/180.0 * M_PI;
      
      //! scaling factor to reduce wind away from center
      double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(xA - cMxA) + SQR(xB - cMxB))) * 1.e-3;
      double wxA = -scale * ReferenceScale::vmax_atm * (cos(alpha) * (xA - cMxA) + sin(alpha) * (xB - cMxB));
      double wxB = -scale * ReferenceScale::vmax_atm *(-sin(alpha) * (xA - cMxA) + cos(alpha) * (xB - cMxB));

      double ex1 = -sin(x); // elon
      double ex2 =  cos(x);

      return wxA * (ex1*PA1 + ex2 * PA2) + wxB * (ex1 * PB1 + ex2 * PB2);

    }
};
struct AtmY : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
      // xA,xB coordinate in Arctic plane
      double X = Nextsim::EarthRadius * cos(x)*cos(y);
      double Y = Nextsim::EarthRadius * sin(x)*cos(y);
      double Z = Nextsim::EarthRadius * sin(y);
      double xA = X * PA1 + Y * PA2 + Z * PA3;
      double xB = X * PB1 + Y * PB2 + Z * PB3;
      

      constexpr double oneday = 24.0 * 60.0 * 60.0;
      //! Center of cyclone (in m) (cmX, cmX
      double cMx = 51200.0 * time/oneday;
      double cMxA = cMx * PA1 + cMx * PA2; // same in the plane-system
      double cMxB = cMx * PA2 + cMx * PA2;

      
      
      double alpha = 72.0/180.0 * M_PI;
      
      //! scaling factor to reduce wind away from center
      double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(xA - cMxA) + SQR(xB - cMxB))) * 1.e-3; 
      double wxA = -scale * ReferenceScale::vmax_atm * (cos(alpha) * (xA - cMxA) + sin(alpha) * (xB - cMxB));
      double wxB = -scale * ReferenceScale::vmax_atm *(-sin(alpha) * (xA - cMxA) + cos(alpha) * (xB - cMxB));
      
      double ey1 = -sin(y)*cos(x);
      double ey2 = -sin(y)*sin(x);
      double ey3 = cos(y);

      return wxA * (ey1*PA1 + ey2 * PA2 + ey3 * PA3) + wxB * (ey1 * PB1 + ey2 * PB2 + ey3 * PB3);
    }
};
class InitialH : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return 0.3;// + 0.005 * (sin(6.e-5 * x) + sin(3.e-5 * y));
    }
};
class InitialA : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const { return 1.0; }
};

template <int CG, int DGadvection>
void run_benchmark()
{
  //! Compose name of output directory and create it
  std::string resultsdir = "BenchmarkOrcaSpherical_" + std::to_string(CG) + "_" + std::to_string(DGadvection);
  std::filesystem::create_directory(resultsdir);

  Nextsim::ParametricMesh smesh(CoordinateSystem);
  smesh.readmesh("../tests/example3-25km_NH.smesh");
  smesh.TransformToRadians();
  smesh.RotatePoleToGreenland();
  
  // output land mask
  Nextsim::DGVector<1> landmask(smesh);
  for (size_t i=0;i<smesh.nelements;++i)
    landmask(i,0) = smesh.landmask[i];
  Nextsim::VTK::write_dg<1>(resultsdir + "/landmask", 0, landmask, smesh);
				   
  
  // output boundary info
  Nextsim::DGVector<1> boundary(smesh);
  for (size_t j=0;j<4;++j)
    for (size_t i=0;i<smesh.dirichlet[j].size();++i)
      boundary(smesh.dirichlet[j][i],0) = 1+j;
  Nextsim::VTK::write_dg<1>(resultsdir + "/boundary", 0, boundary, smesh);
  
  
  
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
  VP.fc = 0.0; // no coriolis!
  
  std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
	    << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta << std::endl
	    << "CG/DG " << CG << "\t" << DGadvection  << std::endl;
  
  //! VTK output
  constexpr double T_vtk = 60.0*60.0;
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
                Nextsim::VTK::write_dg(resultsdir + "/Delta", printstep, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);

                Nextsim::VTK::write_dg(resultsdir + "/H", printstep, H, smesh);
		Nextsim::VTK::write_dg(resultsdir + "/A", printstep, A, smesh);
                Nextsim::VTK::write_dg(resultsdir + "/Vx", printstep, dgtransport.GetVx(), smesh);
		Nextsim::VTK::write_dg(resultsdir + "/Vy", printstep, dgtransport.GetVy(), smesh);
		
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}

int main()
{
  run_benchmark<2, 1>();
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
