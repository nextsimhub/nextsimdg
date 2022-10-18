/*!
 * @file benchmark_mehlmann_meb.cpp
 * @date 15 August 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "ParametricTransport.hpp"

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
constexpr double T = 2.0 * 24 * 60. * 60.; //!< Time horizon 2 days
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
        double cM = 256000. + 51200. * time / oneday;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cM) + SQR(y - cM))) * 1.e-3;

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * ReferenceScale::vmax_atm * (cos(alpha) * (x - cM) + sin(alpha) * (y - cM));
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
        double cM = 256000. + 51200. * time / oneday;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cM) + SQR(y - cM))) * 1.e-3;

        double alpha = 72.0 / 180.0 * M_PI;

        return -scale * ReferenceScale::vmax_atm * (-sin(alpha) * (x - cM) + cos(alpha) * (y - cM));
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

class InitialD : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const { return 0.0; }
};

//////////////////////////////////////////////////

template <int CG, int DGadvection, int DGstress>
void run_benchmark(const std::string meshfile)
{
    //! Define the spatial mesh
    Nextsim::ParametricMesh smesh;
    smesh.readmesh(meshfile);
    size_t NX = smesh.nx;

    //! Compose name of output directory and create it
    std::string resultsdir = "Benchmark_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "_" + std::to_string(DGstress)
        + "__" + std::to_string(NX);
    std::filesystem::create_directory(resultsdir);

    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG, DGstress> momentum(smesh);

    //! define the time mesh
    constexpr double dt_adv = 60.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    constexpr size_t NT_meb = 100; //!< Momentum substeps
    constexpr double dt_mom = dt_adv / NT_meb; //!< Time step of momentum problem

    //! Rheology-Parameters
    Nextsim::MEBParameters Params;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "Momentum subcycling " << NT_meb << "\t dt_momentum " << dt_mom
              << std::endl;

    //! VTK output
    constexpr double T_vtk = 1. * 60.0 * 60.0; // every 4 hours
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
    Nextsim::DGVector<DGadvection> H(smesh), A(smesh), D(smesh); //!< ice height and concentration
    Nextsim::Interpolations::Function2DG(smesh, H, InitialH());
    Nextsim::Interpolations::Function2DG(smesh, A, InitialA());
    Nextsim::Interpolations::Function2DG(smesh, D, InitialD());

    ////////////////////////////////////////////////// i/o of initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh);
    Nextsim::VTK::write_dg(resultsdir + "/A", 0, A, smesh);
    Nextsim::VTK::write_dg(resultsdir + "/H", 0, H, smesh);
    Nextsim::VTK::write_dg(resultsdir + "/D", 0, D, smesh);
    //Nextsim::VTK::write_dg(resultsdir + "/S11", 0, momentum.GetS11(), smesh);
    //Nextsim::VTK::write_dg(resultsdir + "/S12", 0, momentum.GetS12(), smesh);
    //Nextsim::VTK::write_dg(resultsdir + "/S22", 0, momentum.GetS22(), smesh);
    Nextsim::VTK::write_dg(resultsdir + "/Div", 0,
        Nextsim::Tools::TensorInvI(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
    Nextsim::VTK::write_dg(resultsdir + "/Shear", 0,
        Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
    Nextsim::VTK::write_dg(resultsdir + "/sigma_n", 0,
        Nextsim::Tools::TensorInvI(smesh, momentum.GetS11(), momentum.GetS12(), momentum.GetS22()), smesh);
    Nextsim::VTK::write_dg(resultsdir + "/tau", 0,
        Nextsim::Tools::TensorInvII(smesh, momentum.GetS11(), momentum.GetS12(), momentum.GetS22()), smesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    ////////////////////////////////////////////////// Initialize transport
    Nextsim::ParametricTransport<DGadvection> dgtransport(smesh);
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
        ///dgtransport.prepareAdvection(momentum.GetVx(), momentum.GetVy());
        dgtransport.prepareAdvection(momentum.GetAvgSubiterVx(), momentum.GetAvgSubiterVy());

        dgtransport.step(dt_adv, A);
        dgtransport.step(dt_adv, H);
        dgtransport.step(dt_adv, D);

        //! Gauss-point limiting
        Nextsim::LimitMax(A, 1.0);
        Nextsim::LimitMin(A, 0.0);
        Nextsim::LimitMin(H, 0.0);
        Nextsim::LimitMax(D, 1.0-1.e-12);
        Nextsim::LimitMin(D, 0.0);
        Nextsim::GlobalTimer.stop("time loop - advection");

        //////////////////////////////////////////////////
        Nextsim::GlobalTimer.start("time loop - meb");
        momentum.prepareIteration(H, A, D);
        // MEB momentum subcycling
        for (size_t mebstep = 0; mebstep < NT_meb; ++mebstep) {
            momentum.MEBStep(Params, NT_meb, dt_adv, H, A, D);
            // <- MPI
        }

        Nextsim::GlobalTimer.stop("time loop - meb");

        //////////////////////////////////////////////////
        if (WRITE_VTK) // Output
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                int printstep = timestep / NT_vtk + 1.e-4;
                Nextsim::GlobalTimer.start("time loop - i/o");
                Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", printstep, momentum.GetVx(), momentum.GetVy(), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/A", printstep, A, smesh);
                Nextsim::VTK::write_dg(resultsdir + "/H", printstep, H, smesh);
                Nextsim::VTK::write_dg(resultsdir + "/D", printstep, D, smesh);
                //Nextsim::VTK::write_dg(resultsdir + "/S11", printstep, momentum.GetS11(), smesh);
                //Nextsim::VTK::write_dg(resultsdir + "/S12", printstep, momentum.GetS12(), smesh);
                //Nextsim::VTK::write_dg(resultsdir + "/S22", printstep, momentum.GetS22(), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/Div", printstep, Nextsim::Tools::TensorInvI(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/sigma_n", printstep,
                    Nextsim::Tools::TensorInvI(smesh, momentum.GetS11(), momentum.GetS12(), momentum.GetS22()), smesh);
                Nextsim::VTK::write_dg(resultsdir + "/tau", printstep,
                    Nextsim::Tools::TensorInvII(smesh, momentum.GetS11(), momentum.GetS12(), momentum.GetS22()), smesh);
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}

int main()
{
ssssssssssssssssssssssssssssssssssssssssssss
    //run_benchmark<2, 3, 8>("../ParametricMesh/rectangle_256x256.smesh");
    run_benchmark<2, 3, 8>("../ParametricMesh/rectangle_128x128.smesh");

    // std::vector<std::string> meshes;
    // meshes.push_back("../ParametricMesh/rectangle_16x16.smesh");
    // meshes.push_back("../ParametricMesh/rectangle_32x32.smesh");
    // meshes.push_back("../ParametricMesh/rectangle_64x64.smesh");
    // meshes.push_back("../ParametricMesh/rectangle_128x128.smesh");
    // meshes.push_back("../ParametricMesh/rectangle_256x256.smesh");
    // meshes.push_back("../ParametricMesh/rectangle_512x512.smesh");

    // for (const auto& it : meshes) {
    //     run_benchmark<1, 1, 3>(it);
    //     run_benchmark<1, 3, 3>(it);
    //     run_benchmark<1, 6, 3>(it);
    //     run_benchmark<2, 1, 8>(it);
    //     run_benchmark<2, 3, 8>(it);
    //     run_benchmark<2, 6, 8>(it);
    // }
}
