/*!
 * @file benchmark_mehlmann_mevp.cpp
 * @date 24 July 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 * @date 14 Nov 2022
 * @author Timothy Spain <timothy.spain@nersc.no>

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
#include "DGModelArray.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelState.hpp"
#include "include/ParametricGrid.hpp"
#include "include/ParaGridIO.hpp"
#include "include/gridNames.hpp"

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

// Adapted to test the reading and writing of DG data using the NextsimDG infrastructure.

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
//////////////////////////////////////////////////

template <int CG, int DGadvection, int DGstress>
void run_benchmark(const std::string meshfile)
{
    // Set the numbers of components from the template parameters
    Nextsim::ModelArray::setNComponents(Nextsim::ModelArray::Type::DG, DGadvection);
    Nextsim::ModelArray::setNComponents(Nextsim::ModelArray::Type::DGSTRESS, DGstress);
    Nextsim::ModelArray::setNComponents(Nextsim::ModelArray::Type::VERTEX, Nextsim::ModelArray::nCoords);

    Nextsim::ParametricGrid gridIn;
    Nextsim::ParaGridIO* readIO = new Nextsim::ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    Nextsim::ModelState state = gridIn.getModelState("nextsimMEVPstart.nc");
    size_t NX = Nextsim::ModelArray::definedDimensions.at(Nextsim::ModelArray::Dimension::X).length;
    // Let's ASSUME it's square
    //! Define the spatial mesh
    Nextsim::ModelArray& coords = state.data.at(Nextsim::coordsName);
    Nextsim::ParametricMesh smesh(NX, NX, state.data.at(Nextsim::coordsName).data().matrix());

    //! Define the spatial mesh
    Nextsim::ParametricMesh smeshOG;
    smeshOG.readmesh(meshfile);

    // Assert some facts about the vertex positions
    assert(smesh.nx == NX);
    assert(smesh.ny == NX);
    assert(smesh.nelements == NX * NX);
    assert(smesh.nnodes == (NX + 1) * (NX + 1));
    assert(smesh.vertices(0, 0) == smesh.vertices(1, 0));
    assert(smesh.vertices(0, 1) != smesh.vertices(1, 1));
    assert(smesh.vertices(0, 0) != smesh.vertices(NX+1, 0));
    assert(smesh.vertices(0, 1) == smesh.vertices(NX+1, 1));

    //! Compose name of output directory and create it
    std::string resultsdir = "BenchmarkNS_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "_" + std::to_string(DGstress)
        + "__" + std::to_string(NX);
    std::filesystem::create_directory(resultsdir);

    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG, DGstress> momentum(smesh);

    //! define the time mesh
    constexpr double dt_adv = 120.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 1500.0;
    constexpr double beta = 1500.0;
    constexpr size_t NT_evp = 100;

    //! Rheology-Parameters
    Nextsim::VPParameters VP;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta << std::endl
              << "CG/DG/DGstress " << CG << "\t" << DGadvection << "\t" << DGstress
              << std::endl;

    //! VTK output
    constexpr double T_vtk = 48. * 60.0 * 60.0; // only at end
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 10.0 * 60.0; // every 30 minute
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    // netCDF output
    std::string diagFile = "mevp.nc";
    size_t ncOutputPeriod = 30; // Every hour

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

    Nextsim::HField mask(Nextsim::ModelArray::Type::H);
    mask = 1.;

    Nextsim::DGField& msH = state.data.at(Nextsim::hiceName);
    Nextsim::DGField& msA = state.data.at(Nextsim::ciceName);

    Nextsim::DGModelArray::ma2dg(msH, H);
    Nextsim::DGModelArray::ma2dg(msA, A);

    state = {{
            { Nextsim::hiceName, msH },
            { Nextsim::ciceName, msA },
            { Nextsim::coordsName, coords },
            { Nextsim::maskName, mask },
    }, {}
    };

    Nextsim::ModelMetadata metadata;
    Nextsim::TimePoint startTime("2000-01-01T00:00:00Z");
    metadata.setTime(startTime);

    readIO->writeDiagnosticTime(state, metadata, diagFile);

    ////////////////////////////////////////////////// i/o of initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    if (0) // write initial?
        if (WRITE_VTK) {
            Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh);
            Nextsim::VTK::write_dg(resultsdir + "/A", 0, A, smesh);
            Nextsim::VTK::write_dg(resultsdir + "/H", 0, H, smesh);
            Nextsim::VTK::write_dg(resultsdir + "/Shear", 0, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
        }
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
        if (!(timestep % ncOutputPeriod)) {
            Nextsim::DGField ncH(Nextsim::ModelArray::Type::DG);
            Nextsim::DGField ncA(Nextsim::ModelArray::Type::DG);

            Nextsim::DGModelArray::dg2ma<DGadvection>(H, ncH);
            Nextsim::DGModelArray::dg2ma<DGadvection>(A, ncA);

            Nextsim::ModelState state = {{
                    { Nextsim::hiceName, ncH },
                    { Nextsim::ciceName, ncA },
            }, {}
            };

            Nextsim::TimePoint nowTime(startTime + Nextsim::Duration(time));
            metadata.setTime(nowTime);
            readIO->writeDiagnosticTime(state, metadata, diagFile);

        }
    }
    Nextsim::GlobalTimer.stop("time loop");

    readIO->close(diagFile);

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}

int main()
{
    run_benchmark<1, 1, 3>("./distortedrectangle_128x128.smesh");

    // std::vector<std::string> meshes;
    // meshes.push_back("../ParametricMesh/distortedrectangle_16x16.smesh");
    // meshes.push_back("../ParametricMesh/distortedrectangle_32x32.smesh");
    // meshes.push_back("../ParametricMesh/distortedrectangle_64x64.smesh");
    // meshes.push_back("../ParametricMesh/distortedrectangle_128x128.smesh");
    // meshes.push_back("../ParametricMesh/distortedrectangle_256x256.smesh");
    // meshes.push_back("../ParametricMesh/distortedrectangle_512x512.smesh");

    // for (const auto& it : meshes) {
    //     run_benchmark<1, 1, 3>(it);
    //     run_benchmark<1, 3, 3>(it);
    //     run_benchmark<1, 6, 3>(it);
    //     run_benchmark<2, 1, 8>(it);
    //     run_benchmark<2, 3, 8>(it);
    //     run_benchmark<2, 6, 8>(it);
    // }
}
