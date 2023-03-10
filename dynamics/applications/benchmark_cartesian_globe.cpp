/*!
 * @file benchmark_cartesian_globe.cpp
 * @date 24 July 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

/*!
 *
 * Ice dynamics test case on a Cartesian projection
 * of the globe.
 * Domain of size [0,40000] * [-10000,10000]
 * Ice on [0,40000] * [-8000,8000] (by landmask)
 * Periodic at x=0 and x = 40000, Dirichlet at y = +/- 10000
 *
 * Wind is perturbation of (10,0)
 * Ocean (1,0) * (1.0 - (x/8000)^2) and zero on land
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

Nextsim::COORDINATES CoordinateSystem = Nextsim::CARTESIAN;
namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return x * x; }

//////////////////////////////////////////////////// Benchmark testcase from [Mehlmann / Richter, ...]
//! Description of the problem data, wind & ocean fields

namespace ReferenceScale {
constexpr double T = 30.0 * 24 * 60. * 60.; //!< Time horizon 1 month

constexpr double Lx = 40000.e3; //!< Size of domain, [0,Lx] x [-Ly,Ly]
constexpr double Ly = 10000.e3;
}

class OceanX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      double r2 = 1.0-pow(y/8000.0e3,2.0);
      if (r2<0.0)
	return 0.0;
      return r2;
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
      return 10.0 + 5.0 * sin(4.0*M_PI*x/ReferenceScale::Lx) * cos(4.0*M_PI*y/ReferenceScale::Ly);
    }
};

struct AtmY : public Nextsim::Interpolations::Function {

public:
  double time;
    double operator()(double x, double y) const
    {
      return 5.0 * sin(4.0*M_PI*(x/ReferenceScale::Lx-time/ReferenceScale::T)) * cos(4.0*M_PI*y/ReferenceScale::Ly);
    }
};
class InitialH : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return 0.5;
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

void create_mesh(const std::string& meshname, const size_t N)
{
    std::ofstream OUT(meshname.c_str());
    OUT << "ParametricMesh 2.0" << std::endl
        << 2*N << "\t" << N << std::endl;
    for (size_t iy = 0; iy <= N; ++iy)
        for (size_t ix = 0; ix <= 2*N; ++ix)
	  OUT <<  ReferenceScale::Lx * ix / (2*N) << "\t"
	      << -ReferenceScale::Ly + 2.0 * ReferenceScale::Ly * iy / N << std::endl;


    OUT << "landmask " << 2*N*N << std::endl;
    for (size_t iy = 0; iy < N; ++iy)
      {
	for (size_t ix = 0; ix < 2*N; ++ix)
	  if ( (iy<0.1*N) || (iy>=0.9*N-1.e-8) )
	    OUT << "0 ";
	  else
	    OUT << "1 ";
	OUT << std::endl;
      }

    // dirichlet top and bottom not required due to land mask
    OUT << "dirichlet 0" << std::endl;

    // periodic left/right
    OUT << "periodic 1" << std::endl;

    OUT << N << std::endl;
    for (size_t i = 0; i < N; ++i)
      OUT << 2*N - 1 + i * 2 * N << "\t" << i * 2 * N << "\t1" << std::endl;

    OUT.close();
}
template <int CG, int DGadvection>
void run_benchmark(const size_t N)
{
    //! Define the spatial mesh
    create_mesh("tmp-benchmark.smesh", N);
    Nextsim::ParametricMesh smesh(CoordinateSystem);
    smesh.readmesh("tmp-benchmark.smesh");

    //! Compose name of output directory and create it
    std::string resultsdir = "CartesianGlobe_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(N);
    std::filesystem::create_directory(resultsdir);

    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh);

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

    ////////////////////////////////////////////////// Forcing
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceanx(), OceanX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceany(), OceanY());
    AtmX AirX;  AirX.time = 0.0;
    AtmY AirY;  AirY.time = 0.0;
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AirX);
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AirY);



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

        if (timestep % NT_log == 0)
            std::cout << "\rAdvection step " << timestep << "\t " << std::setprecision(2)
                      << std::fixed << std::setw(10) << std::right << time << "s\t" << std::setw(8)
                      << std::right << timeInMinutes << "m\t" << std::setw(6) << std::right
                      << timeInHours << "h\t" << std::setw(6) << std::right << timeInDays << "d\t\t"
                      << std::flush;

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
  run_benchmark<2, 6>(64);
}
