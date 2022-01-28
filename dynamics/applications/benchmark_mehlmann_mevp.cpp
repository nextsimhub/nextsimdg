#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cgmomentum.hpp"
#include "cgvector.hpp"
#include "dginitial.hpp"
#include "dgvisu.hpp"
#include "dynamics.hpp"
#include "mevp.hpp"
#include "stopwatch.hpp"
#include "tools.hpp"

bool WRITE_VTK = true;

#define CG 1
#define DGadvection 0
#define DGstress 0

namespace Nextsim {
extern Timer GlobalTimer;
}

namespace ReferenceScale {
// Benchmark testcase from [Mehlmann / Richter, ...]
constexpr double T = 2 * 24. * 60. * 60.; //!< Time horizon 2 days
constexpr double L = 512000.0; //!< Size of domain
constexpr double vmax_ocean = 0.01; //!< Maximum velocity of ocean
constexpr double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind

constexpr double rho_ice = 900.0; //!< Sea ice density
constexpr double rho_atm = 1.3; //!< Air density
constexpr double rho_ocean = 1026.0; //!< Ocean density

constexpr double C_atm = 1.2e-3; //!< Air drag coefficient
constexpr double C_ocean = 5.5e-3; //!< Ocean drag coefficient

constexpr double F_atm = C_atm * rho_atm; //!< effective factor for atm-forcing
constexpr double F_ocean = C_ocean * rho_ocean; //!< effective factor for ocean-forcing

constexpr double Pstar = 27500.0; //!< Ice strength
constexpr double fc = 1.46e-4; //!< Coriolis

constexpr double DeltaMin = 2.e-9; //!< Viscous regime
}

inline constexpr double SQR(double x)
{
    return x * x;
}

//! Description of the problem data, wind & ocean fields
class OceanX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (2.0 * y / ReferenceScale::L - 1.0);
    }
};
class OceanY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (1.0 - 2.0 * x / ReferenceScale::L);
    }
};

class AtmX : virtual public Nextsim::InitialBase {
    double time;

public:
    void settime(double t)
    {
        time = t;
    }
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
class AtmY : virtual public Nextsim::InitialBase {
    double time;

public:
    void settime(double t)
    {
        time = t;
    }
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
class InitialH : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 0.3 + 0.005 * (sin(6.e-5 * x) + sin(3.e-5 * y));
    }
};
class InitialA : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 1.0;
    }
};

int main()
{
    //!

    Nextsim::CGMomentum momentum;

    //! Define the spatial mesh
    Nextsim::Mesh mesh;
    constexpr size_t N = 128; //!< Number of mesh nodes
    mesh.BasicInit(N, N, ReferenceScale::L / N, ReferenceScale::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    //! define the time mesh
    constexpr double dt_adv = 120.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 800.0;
    constexpr double beta = 800.0;
    constexpr size_t NT_evp = 100;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta
              << std::endl;

    //! VTK output
    constexpr double T_vtk = 1.0 * 60.0 * 60.0; // evey 4 hours
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 10.0 * 60.0; // every 30 minute
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    //! Variables
    Nextsim::CGVector<CG> vx(mesh), vy(mesh); //!< velocity
    Nextsim::CGVector<CG> tmpx(mesh), tmpy(mesh); //!< tmp for stress.. should be removed
    Nextsim::CGVector<CG> cg_A(mesh), cg_H(mesh); //!< interpolation of ice height and conc.
    Nextsim::CGVector<CG> vx_mevp(mesh), vy_mevp(mesh); //!< temp. Velocity used for MEVP
    vx.zero();
    vy.zero();
    vx_mevp = vx;
    vy_mevp = vy;

    Nextsim::CGVector<CG> OX(mesh); //!< x-component of ocean velocity
    Nextsim::CGVector<CG> OY(mesh); //!< y-component of ocean velocity
    Nextsim::InterpolateCG(mesh, OX, OceanX());
    Nextsim::InterpolateCG(mesh, OY, OceanY());
    Nextsim::CGVector<CG> AX(mesh); //!< x-component of atm. velocity
    Nextsim::CGVector<CG> AY(mesh); //!< y-component of atm. velocity
    AtmX AtmForcingX;
    AtmY AtmForcingY;
    AtmForcingX.settime(0.0);
    AtmForcingY.settime(0.0);
    AX.zero();
    AY.zero();
    // Nextsim::InterpolateCG(mesh, AX, AtmForcingX);
    // Nextsim::InterpolateCG(mesh, AY, AtmForcingY);

    Nextsim::CellVector<DGadvection> H(mesh), A(mesh); //!< ice height and concentration
    Nextsim::L2ProjectInitial(mesh, H, InitialH());
    Nextsim::L2ProjectInitial(mesh, A, InitialA());
    Nextsim::CellVector<DGstress> E11(mesh), E12(mesh), E22(mesh); //!< storing strain rates
    Nextsim::CellVector<DGstress> S11(mesh), S12(mesh), S22(mesh); //!< storing stresses rates

    Nextsim::CellVector<0> DELTA(mesh); //!< Storing DELTA
    Nextsim::CellVector<0> SHEAR(mesh); //!< Storing DELTA
    Nextsim::CellVector<0> S1(mesh), S2(mesh); //!< Stress invariants

    // save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg("ResultsBenchmark/vx", 0, vx, mesh);
    Nextsim::VTK::write_cg("ResultsBenchmark/vy", 0, vy, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/A", 0, A, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/H", 0, H, mesh);

    Nextsim::Tools::Delta(mesh, E11, E12, E22, ReferenceScale::DeltaMin, DELTA);
    Nextsim::VTK::write_dg("ResultsBenchmark/Delta", 0, DELTA, mesh);
    Nextsim::Tools::Shear(mesh, E11, E12, E22, ReferenceScale::DeltaMin, SHEAR);
    Nextsim::VTK::write_dg("ResultsBenchmark/Shear", 0, SHEAR, mesh);

    Nextsim::VTK::write_dg("ResultsBenchmark/S11", 0, S11, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/S12", 0, S12, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/S22", 0, S22, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/E11", 0, E11, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/E12", 0, E12, mesh);
    Nextsim::VTK::write_dg("ResultsBenchmark/E22", 0, E22, mesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    //! Transport
    Nextsim::CellVector<DGadvection> dgvx(mesh), dgvy(mesh);
    Nextsim::DGTransport<DGadvection> dgtransport(dgvx, dgvy);
    dgtransport.settimesteppingscheme("rk2");
    dgtransport.setmesh(mesh);

    //! Initial Forcing
    AtmForcingX.settime(0);
    AtmForcingY.settime(0);
    Nextsim::InterpolateCG(mesh, AX, AtmForcingX);
    Nextsim::InterpolateCG(mesh, AY, AtmForcingY);

    Nextsim::GlobalTimer.start("time loop");

    for (size_t timestep = 1; timestep <= NT; ++timestep) {
        double time = dt_adv * timestep;
        double timeInMinutes = time / 60.0;
        double timeInHours = time / 60.0 / 60.0;
        double timeInDays = time / 60.0 / 60.0 / 24.;

        if (timestep % NT_log == 0)
            std::cout << "\rAdvection step " << timestep << "\t "
                      << std::setprecision(2)
                      << std::fixed
                      << std::setw(10) << std::right
                      << time << "s\t"
                      << std::setw(8) << std::right
                      << timeInMinutes << "m\t"
                      << std::setw(6) << std::right
                      << timeInHours << "h\t"
                      << std::setw(6) << std::right
                      << timeInDays << "d\t\t" << std::flush;

        //! Initialize time-dependent data
        Nextsim::GlobalTimer.start("time loop - forcing");
        AtmForcingX.settime(time);
        AtmForcingY.settime(time);
        Nextsim::InterpolateCG(mesh, AX, AtmForcingX);
        Nextsim::InterpolateCG(mesh, AY, AtmForcingY);
        Nextsim::GlobalTimer.stop("time loop - forcing");

        //! Advection step
        Nextsim::GlobalTimer.start("time loop - advection");
        momentum.ProjectCGToDG(mesh, dgvx, vx);
        momentum.ProjectCGToDG(mesh, dgvy, vy);
        dgtransport.reinitvelocity();
        dgtransport.step(dt_adv, A);
        dgtransport.step(dt_adv, H);

        //! Very simple limiting (just constants)
        A.col(0) = A.col(0).cwiseMin(1.0);
        A.col(0) = A.col(0).cwiseMax(0.0);
        H.col(0) = H.col(0).cwiseMax(0.0);

        momentum.InterpolateDGToCG(mesh, cg_A, A);
        momentum.InterpolateDGToCG(mesh, cg_H, H);

        cg_A = cg_A.cwiseMin(1.0);
        cg_A = cg_A.cwiseMax(0.0);
        cg_H = cg_H.cwiseMax(1.e-4); //!< Limit H from below

        Nextsim::GlobalTimer.stop("time loop - advection");

        Nextsim::GlobalTimer.start("time loop - mevp");
        //! Store last velocity for MEVP
        vx_mevp = vx;
        vy_mevp = vy;

        //! MEVP subcycling
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {

            Nextsim::GlobalTimer.start("time loop - mevp - strain");
            //! Compute Strain Rate
            momentum.ProjectCG2VelocityToDG1Strain(mesh, E11, E12, E22, vx, vy);
            Nextsim::GlobalTimer.stop("time loop - mevp - strain");

            Nextsim::GlobalTimer.start("time loop - mevp - stress");

            Nextsim::mEVP::StressUpdate(mesh, S11, S12, S22,
                E11, E12, E22, H, A,
                ReferenceScale::Pstar,
                ReferenceScale::DeltaMin,
                alpha, beta);

            Nextsim::GlobalTimer.stop("time loop - mevp - stress");

            Nextsim::GlobalTimer.start("time loop - mevp - update");
            //! Update
            Nextsim::GlobalTimer.start("time loop - mevp - update1");

            //	    update by a loop.. implicit parts and h-dependent

            //
            vx = (1.0 / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                      + cg_A.array() * ReferenceScale::F_ocean * (OX.array() - vx.array()).abs()) // implicit parts
                * (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (beta * vx.array() + vx_mevp.array()) + // pseudo-timestepping
                    cg_A.array() * (ReferenceScale::F_atm * AX.array().abs() * AX.array() + // atm forcing
                        ReferenceScale::F_ocean * (OX - vx).array().abs() * OX.array()) // ocean forcing
                    + ReferenceScale::rho_ice * cg_H.array() * ReferenceScale::fc * (vy - OY).array() // cor + surface
                    ))
                     .matrix();
            vy = (1.0 / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                      + cg_A.array() * ReferenceScale::F_ocean * (OY.array() - vy.array()).abs()) // implicit parts
                * (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (beta * vy.array() + vy_mevp.array()) + // pseudo-timestepping
                    cg_A.array() * (ReferenceScale::F_atm * AY.array().abs() * AY.array() + // atm forcing
                        ReferenceScale::F_ocean * (OY - vy).array().abs() * OY.array()) // ocean forcing
                    + ReferenceScale::rho_ice * cg_H.array() * ReferenceScale::fc * (OX - vx).array() // cor + surface
                    ))
                     .matrix();
            Nextsim::GlobalTimer.stop("time loop - mevp - update1");

            Nextsim::GlobalTimer.start("time loop - mevp - update2");
            // Implicit etwas ineffizient
            tmpx.zero();
            tmpy.zero();
            momentum.AddStressTensor(mesh, -1.0, tmpx, tmpy, S11, S12, S22);
            vx += (1.0 / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                       + cg_A.array() * ReferenceScale::F_ocean * (OX.array() - vx.array()).abs()) // implicit parts
                * tmpx.array())
                      .matrix();
            vy += (1.0 / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                       + cg_A.array() * ReferenceScale::F_ocean * (OY.array() - vy.array()).abs()) // implicit parts
                * tmpy.array())
                      .matrix();

            Nextsim::GlobalTimer.stop("time loop - mevp - update2");
            Nextsim::GlobalTimer.stop("time loop - mevp - update");

            Nextsim::GlobalTimer.start("time loop - mevp - bound.");
            momentum.DirichletZero(mesh, vx);
            momentum.DirichletZero(mesh, vy);
            Nextsim::GlobalTimer.stop("time loop - mevp - bound.");
        }
        Nextsim::GlobalTimer.stop("time loop - mevp");

        //         //! Output
        if (WRITE_VTK)
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                size_t printstep = timestep / NT_vtk + 1.e-4;

                char s[80];
                sprintf(s, "ResultsBenchmark/ellipse_%03d.txt");
                std::ofstream ELLOUT(s);
                for (size_t i = 0; i < mesh.n; ++i)
                    ELLOUT << S1(i, 0) << "\t" << S2(i, 0) << std::endl;
                ELLOUT.close();

                Nextsim::GlobalTimer.start("time loop - i/o");
                Nextsim::VTK::write_cg("ResultsBenchmark/vx", printstep, vx, mesh);
                Nextsim::VTK::write_cg("ResultsBenchmark/vy", printstep, vy, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/A", printstep, A, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/H", printstep, H, mesh);

                Nextsim::Tools::Delta(mesh, E11, E12, E22, ReferenceScale::DeltaMin, DELTA);
                Nextsim::VTK::write_dg("ResultsBenchmark/Delta", printstep, DELTA, mesh);
                Nextsim::Tools::Shear(mesh, E11, E12, E22, ReferenceScale::DeltaMin, SHEAR);
                Nextsim::VTK::write_dg("ResultsBenchmark/Shear", printstep, SHEAR, mesh);

                Nextsim::VTK::write_dg("ResultsBenchmark/S11", printstep, S11, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/S12", printstep, S12, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/S22", printstep, S22, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/E11", printstep, E11, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/E12", printstep, E12, mesh);
                Nextsim::VTK::write_dg("ResultsBenchmark/E22", printstep, E22, mesh);

                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
