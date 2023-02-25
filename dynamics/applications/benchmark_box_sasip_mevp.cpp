/*!
 * @file benchmark_mehlmann_mevp.cpp
 * @date 24 July 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "ParametricTransport.hpp"
#include "Tools.hpp"
#include "VPParameters.hpp"
#include "cgParametricMomentum.hpp"
#include "cgVector.hpp"
#include "dgInitial.hpp"
#include "dgLimit.hpp"
#include "dgVisu.hpp"

#include "stopwatch.hpp"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

bool WRITE_VTK = true;

Nextsim::COORDINATES CoordinateSystem = Nextsim::CARTESIAN;
/*!
 *
 * Sets the order of the velocity (CG) of advection (DGadvection) and
 * of the stress & strain. This should give the gradient space of the
 * CG space for stability. CG=1 -> DGstress=3, CG=2 -> DGstress -> 8
 */
#define CG 2
#define DGadvection 3

namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return x * x; }

//////////////////////////////////////////////////// Benchmark testcase from [Mehlmann / Richter, ...]
//! Description of the problem data, wind & ocean fields

namespace ReferenceScale {
constexpr double T = 30. * 24 * 60. * 60.; //!< Time horizon 2 days
constexpr double L = 1000000.0; //!< Size of domain !!!
constexpr double vmax_ocean = 0.1; //!< Maximum velocity of ocean
const double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind
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

class AtmX : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
        double X = M_PI * x / ReferenceScale::L;
        double Y = M_PI * y / ReferenceScale::L;
        constexpr double T = 4.0 * 24.0 * 60.0 * 60.0; //!< 4 days
        return 5.0 + (sin(2 * M_PI * time / T) - 3.0) * sin(2 * X) * sin(Y);
    }
};
class AtmY : public Nextsim::Interpolations::Function {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
        double X = M_PI * x / ReferenceScale::L;
        double Y = M_PI * y / ReferenceScale::L;
        constexpr double T = 4.0 * 24.0 * 60.0 * 60.0; //!< 4 days
        return 5.0 + (sin(2 * M_PI * time / T) - 3) * sin(2 * Y) * sin(X);
    }
};

class InitialH : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return 2.0; //!< constant ice height for box test
    }
};
class InitialA : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const { return x / ReferenceScale::L; }
};


//////////////////////////////////////////////////
void create_mesh(const std::string& meshname, const size_t Nx, const double distort)
{
    std::ofstream OUT(meshname.c_str());
    OUT << "ParametricMesh 2.0" << std::endl
        << Nx << "\t" << Nx << std::endl;
    for (size_t iy = 0; iy <= Nx; ++iy)
        for (size_t ix = 0; ix <= Nx; ++ix)
            OUT << ReferenceScale::L * ix / Nx + ReferenceScale::L * distort * sin(M_PI * ix / Nx * 3.0) * sin(M_PI * iy / Nx) << "\t"
                << ReferenceScale::L * iy / Nx + ReferenceScale::L * distort * sin(M_PI * iy / Nx * 2.0) * sin(M_PI * ix / Nx * 2.0) << std::endl;


    OUT << "landmask 0" << std::endl;
    
    if (1) // std-setup, dirichlet everywhere
    {
        // dirichlet info
      OUT << "dirichlet " << 2*Nx+2*Nx << std::endl;
      for (size_t i = 0; i < Nx; ++i)
	OUT << i << "\t" << 0 << std::endl; // lower
      for (size_t i = 0; i < Nx; ++i)
	OUT << Nx * (Nx - 1) + i << "\t" << 2 << std::endl; // upper
      
      for (size_t i = 0; i < Nx; ++i)
	OUT << i * Nx << "\t" << 3 << std::endl; // left
      for (size_t i = 0; i < Nx; ++i)
	OUT << i * Nx + Nx - 1 << "\t" << 1 << std::endl; // right
      
        OUT << "periodic 0" << std::endl;
    } else {
        OUT << "dirichlet 0" << std::endl;
        OUT << "periodic 2" << std::endl;

        OUT << 2 * Nx << std::endl; // horizontal (x-)
        for (size_t i = 0; i < Nx; ++i)
            OUT << Nx * (Nx - 1) + i << "\t" << i << "\t0" << std::endl; // left
        for (size_t i = 0; i < Nx; ++i)
            OUT << Nx - 1 + i * Nx << "\t" << i * Nx << "\t1" << std::endl;
    }
    OUT.close();
}

int main()
{
  const int NX = 32;
    //! Define the spatial mesh
  create_mesh("benchmark_box.smesh",NX,0.0);

  std::string resultsdir = "BenchmarkBox_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(NX);
  std::filesystem::create_directory(resultsdir);

  Nextsim::ParametricMesh smesh(CoordinateSystem);
    smesh.readmesh("benchmark_box.smesh");

    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh);

    //! define the time mesh
    constexpr double dt_adv = 120.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 600.0;
    constexpr double beta = 600.0;
    constexpr size_t NT_evp = 200;
    Nextsim::VPParameters VP;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta
              << std::endl;

    //! VTK output
    constexpr double T_vtk = 1.0 * 24.0 * 60.0 * 60.0; // every day
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 10.0 * 60.0; // every 10 minute
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    ////////////////////////////////////////////////// Forcing
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceanx(), OceanX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetOceany(), OceanY());

    AtmX AtmForcingX;
    AtmY AtmForcingY;
    AtmForcingX.settime(0.0);
    AtmForcingY.settime(0.0);

    ////////////////////////////////////////////////// Variables and Initial Values
    Nextsim::DGVector<DGadvection> H(smesh), A(smesh); //!< ice height and concentration
    Nextsim::Interpolations::Function2DG(smesh, H, InitialH());
    Nextsim::Interpolations::Function2DG(smesh, A, InitialA());

    // i/o of initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg_velocity(resultsdir+"/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh);
    Nextsim::VTK::write_dg(resultsdir+"/A", 0, A, smesh);
    Nextsim::VTK::write_dg(resultsdir+"/H", 0, H, smesh);
    Nextsim::VTK::write_dg(resultsdir+"/Delta", 0, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh);
    Nextsim::VTK::write_dg(resultsdir+"/Shear", 0, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    ////////////////////////////////////////////////// Initialize transport
    Nextsim::ParametricTransport<DGadvection> dgtransport(smesh);
    dgtransport.settimesteppingscheme("rk2");

    //! Initial Forcing
    AtmForcingX.settime(0);
    AtmForcingY.settime(0);
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmx(), AtmForcingX);
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetAtmy(), AtmForcingY);

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
        Nextsim::Interpolations::CG2DG(smesh, dgtransport.GetVx(), momentum.GetVx());
        Nextsim::Interpolations::CG2DG(smesh, dgtransport.GetVy(), momentum.GetVy());

        dgtransport.reinitnormalvelocity();
        dgtransport.step(dt_adv, A);
        dgtransport.step(dt_adv, H);

        //! Gauss-point limiting
        Nextsim::LimitMax(A, 1.0);
        Nextsim::LimitMin(A, 0.0);
        Nextsim::LimitMin(H, 0.0);

        Nextsim::GlobalTimer.stop("time loop - advection");
        //////////////////////////////////////////////////

        //////////////////////////////////////////////////
        Nextsim::GlobalTimer.start("time loop - mevp");
        momentum.prepareIteration(H, A);
        // MEVP subcycling
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {
            momentum.mEVPStep(VP, NT_evp, alpha, beta, dt_adv, H, A);
            // <- MPI
        }
        Nextsim::GlobalTimer.stop("time loop - mevp");

        //////////////////////////////////////////////////
        if (WRITE_VTK) // Output
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                int printstep = timestep / NT_vtk + 1.e-4;
                Nextsim::GlobalTimer.start("time loop - i/o");
                Nextsim::VTK::write_cg_velocity(resultsdir+"/vel", printstep, momentum.GetVx(), momentum.GetVy(), smesh);
                Nextsim::VTK::write_dg(resultsdir+"/A", printstep, A, smesh);
                Nextsim::VTK::write_dg(resultsdir+"/H", printstep, H, smesh);
                Nextsim::VTK::write_dg(resultsdir+"/Delta", printstep, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh);
                Nextsim::VTK::write_dg(resultsdir+"/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
