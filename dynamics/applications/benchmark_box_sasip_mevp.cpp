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
#include "dgTransport.hpp"
#include "dgVisu.hpp"

#include "stopwatch.hpp"

#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

bool WRITE_VTK = true;

/*!
 *
 * Sets the order of the velocity (CG) of advection (DGadvection) and
 * of the stress & strain. This should give the gradient space of the
 * CG space for stability. CG=1 -> DGstress=3, CG=2 -> DGstress -> 8
 */
#define CG 2
#define DGadvection 3
#define DGstress 8
#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))

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

int main()
{
    //! Define the spatial mesh
    Nextsim::ParametricMesh smesh;
    smesh.readmesh("../ParametricMesh/distortedbox.smesh");

    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG, DGstress> momentum(smesh);

    //! define the time mesh
    constexpr double dt_adv = 120.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 300.0;
    constexpr double beta = 300.0;
    constexpr size_t NT_evp = 100;
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
    Nextsim::CellVector<DGadvection> H(smesh), A(smesh); //!< ice height and concentration
    Nextsim::Interpolations::Function2DG(smesh, H, InitialH());
    Nextsim::Interpolations::Function2DG(smesh, A, InitialA());

    // i/o of initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg_velocity("ResultsBoxParametricMesh/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh);
    Nextsim::VTK::write_dg("ResultsBoxParametricMesh/A", 0, A, smesh);
    Nextsim::VTK::write_dg("ResultsBoxParametricMesh/H", 0, H, smesh);
    Nextsim::VTK::write_dg("ResultsBoxParametricMesh/Delta", 0, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh);
    Nextsim::VTK::write_dg("ResultsBoxParametricMesh/Shear", 0, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    ////////////////////////////////////////////////// Initialize transport
    Nextsim::ParametricTransport<DGadvection, EDGEDOFS(DGadvection)> dgtransport(smesh);
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
        momentum.mEVPIteration(VP, NT_evp, alpha, beta, dt_adv, H, A);
        Nextsim::GlobalTimer.stop("time loop - mevp");

        //////////////////////////////////////////////////
        if (WRITE_VTK) // Output
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                int printstep = timestep / NT_vtk + 1.e-4;
                Nextsim::GlobalTimer.start("time loop - i/o");
                Nextsim::VTK::write_cg_velocity("ResultsBoxParametricMesh/vel", printstep, momentum.GetVx(), momentum.GetVy(), smesh);
                Nextsim::VTK::write_dg("ResultsBoxParametricMesh/A", printstep, A, smesh);
                Nextsim::VTK::write_dg("ResultsBoxParametricMesh/H", printstep, H, smesh);
                Nextsim::VTK::write_dg("ResultsBoxParametricMesh/Delta", printstep, Nextsim::Tools::Delta(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22(), VP.DeltaMin), smesh);
                Nextsim::VTK::write_dg("ResultsBoxParametricMesh/Shear", printstep, Nextsim::Tools::Shear(smesh, momentum.GetE11(), momentum.GetE12(), momentum.GetE22()), smesh);
                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
