/*!
 * @file benchmark_MEB_CGvel.cpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.no>
 */

#include "CheckPoints.hpp"
#include "MEB.hpp"
#include "MEBSandbox.hpp"
#include "Tools.hpp"
#include "cgMomentum.hpp"
#include "cgVector.hpp"
#include "dgInitial.hpp"
#include "dgLimit.hpp"
#include "dgTransport.hpp"
#include "dgVisu.hpp"
#include "mevp.hpp"
#include "stopwatch.hpp"

#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

bool WRITE_VTK = true;

#define CG 1
#define DGadvection 1
#define DGstress 1

#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))

namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return x * x; }

//! Description of the problem data, wind & ocean fields
struct OceanX {
public:
    double operator()(double x, double y) const
    {
        //return RefScaleCanada::vmax_ocean * (2.0 * y / RefScaleCanada::L - 1.0);
        return 0.0;
    }
};
struct OceanY {
public:
    double operator()(double x, double y) const
    {
        //return RefScaleCanada::vmax_ocean * (1.0 - 2.0 * x / RefScaleCanada::L);
        return 0.0;
    }
};

struct AtmX {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
        return 0.0;
    }
};
struct AtmY {
    double time;

public:
    void settime(double t) { time = t; }
    double operator()(double x, double y) const
    {
        constexpr double twohours = 2.0 * 60.0 * 60.0;
        //! maximum is 20 m/s corresponds to 0.625 Pa forcing
        if (time < twohours)
            return -20. * 0.5 * (1. - cos(M_PI * time / twohours));
        else
            return -20.;
    }
};
struct InitialH {
public:
    double operator()(double x, double y) const { return 1.0; }
};
struct InitialA {
public:
    double operator()(double x, double y) const { return 1.0; }
};
struct InitialD {
public:
    double operator()(double x, double y) const { return 0.0; }
};

int main()
{
    //!
    Nextsim::CGMomentum momentum;

    //! Define the spatial mesh
    Nextsim::Mesh mesh;
    constexpr size_t N = 60; //!< Number of mesh nodes
    //restricted domain for the sake of tests
    mesh.BasicInit(N, 2 * N, RefScaleCanada::L / N, RefScaleCanada::L / N);
    //mesh.BasicInit(N, 250, RefScaleCanada::L / N, RefScaleCanada::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    //! define the time mesh
    constexpr double dt_adv = .01; //!< Time step of advection problem
    constexpr size_t NT = RefScaleCanada::T / dt_adv + 1.e-4; //!< Number of Advections steps

    constexpr size_t mom_substeps = 100;
    constexpr double dt_momentum = dt_adv / mom_substeps; //!< Time step of momentum problem

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "Momentum subcycling " << mom_substeps << "\t dt_momentum " << dt_momentum
              << std::endl;

    //! VTK output
    constexpr double T_vtk = .1 * 60.0 * 60.0; // evey 4 hours
    constexpr size_t NT_vtk = T_vtk / dt_adv + 1.e-4;
    //! LOG message
    constexpr double T_log = 10.0 * 60.0; // every 30 minute
    constexpr size_t NT_log = T_log / dt_adv + 1.e-4;

    //! Variables
    Nextsim::CGVector<CG> vx(mesh), vy(mesh); //!< velocity
    Nextsim::CGVector<CG> tmpx(mesh), tmpy(mesh); //!< tmp for stress.. should be removed
    Nextsim::CGVector<CG> cg_A(mesh), cg_H(mesh); //!< interpolation of ice height and conc.
    // Nextsim::CGVector<CG> vx_mevp(mesh), vy_mevp(mesh); //!< temp. Velocity used for MEVP
    vx.zero();
    vy.zero();

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

    Nextsim::CellVector<DGadvection> H(mesh), A(mesh); //!< ice height and concentration
    Nextsim::L2ProjectInitial(mesh, H, InitialH());
    Nextsim::L2ProjectInitial(mesh, A, InitialA());
    Nextsim::CellVector<DGstress> E11(mesh), E12(mesh), E22(mesh); //!< storing strain rates
    Nextsim::CellVector<DGstress> S11(mesh), S12(mesh), S22(mesh); //!< storing stresses rates

    Nextsim::CellVector<1> DELTA(mesh); //!< Storing DELTA
    Nextsim::CellVector<1> SHEAR(mesh); //!< Storing DELTA
    Nextsim::CellVector<DGstress> Sinv1(mesh), Sinv2(mesh); //!< Stress invariants
    Nextsim::CellVector<DGstress> Einv1(mesh), Einv2(mesh); //!< Strain invariants

    // Temporary variables
    Nextsim::CellVector<1> eta1(mesh), eta2(mesh);
    Nextsim::CellVector<1> MU1(mesh), MU2(mesh);
    Nextsim::CellVector<1> stressrelax(mesh);
    Nextsim::CellVector<1> sigma_outside(mesh);
    Nextsim::CellVector<1> tildeP(mesh);
    Nextsim::CellVector<1> Pmax(mesh);
    Nextsim::CellVector<1> td(mesh);

    Nextsim::CellVector<1> tau(mesh);
    Nextsim::CellVector<1> sigma_n(mesh);
    Nextsim::CellVector<1> Regime(mesh);
    Nextsim::CellVector<1> Multip(mesh);
    Nextsim::CellVector<1> Lambda(mesh);

    Nextsim::CellVector<DGadvection> d_crit(mesh);
    Nextsim::CellVector<DGadvection> D(mesh); //!< ice damage. ?? Really dG(0) ??
    Nextsim::L2ProjectInitial(mesh, D, InitialD());

    //! Transport
    Nextsim::CellVector<DGadvection> dgvx(mesh), dgvy(mesh);
    Nextsim::DGTransport<DGadvection, EDGEDOFS(DGadvection)> dgtransport(dgvx, dgvy);
    dgtransport.settimesteppingscheme("rk1");
    dgtransport.setmesh(mesh);

    //read initial
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/vy.00014.txt", vy);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/vx.00014.txt", vx);
    ////Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/A.00014.txt", A);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/H.00014.txt", H);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/D.00014.txt", D);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/S12.00014.txt", S12);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/S11.00014.txt", S11);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/S22.00014.txt", S22);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/E11.00014.txt", E11);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/E12.00014.txt", E12);
    //Nextsim::CheckPoints::loadData("RestartCompression/Checkpoints/E22.00014.txt", E22);

    // save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg("ResultsCompression/vx", 0, vx, mesh);
    Nextsim::VTK::write_cg("ResultsCompression/vy", 0, vy, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/A", 0, A, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/H", 0, H, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/D", 0, D, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/S11", 0, S11, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/S12", 0, S12, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/S22", 0, S22, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/dvx", 0, dgvx, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/dvy", 0, dgvy, mesh);

    Nextsim::Tools::TensorInvariants(mesh, E11, E12, E22, Einv1, Einv2);
    Nextsim::VTK::write_dg("ResultsCompression/Einv1", 0, Einv1, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/Einv2", 0, Einv2, mesh);

    Nextsim::Tools::TensorInvariants(mesh, S11, S12, S22, Sinv1, Sinv2);
    Nextsim::VTK::write_dg("ResultsCompression/sigma_n", 0, Sinv1, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/tau", 0, Sinv2, mesh);

    //Nextsim::Tools::ElastoParams(
    //    mesh, E11, E12, E22, H, A, RefScaleCanada::DeltaMin, RefScaleCanada::Pstar, MU1, MU2);
    //
    //Nextsim::VTK::write_dg("ResultsCompression/mu1", 0, MU1, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/mu2", 0, MU2, mesh);

    //Nextsim::VTK::write_dg("ResultsCompression/eta1", 0, eta1, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/eta2", 0, eta2, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/stressrelax", 0, stressrelax, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/stressOutside", 0, sigma_outside, mesh);
    //Nextsim::VTK::write_dg("ResultsCompression/Pmax", 0, Pmax, mesh);
    //Nextsim::VTK::write_dg("ResultsMEB_CGvel/sigma_n", 0, sigma_n, mesh);
    //Nextsim::VTK::write_dg("ResultsMEB_CGvel/tau", 0, tau, mesh);
    //Nextsim::VTK::write_dg("ResultsMEB_CGvel/td", 0, td, mesh);
    Nextsim::VTK::write_dg("ResultsCompression/d_crit", 0, d_crit, mesh);
    //Nextsim::VTK::write_dg("ResultsMEB_CGvel/Regime", 0, Regime, mesh);
    //Nextsim::VTK::write_dg("ResultsMEB_CGvel/Multip", 0, Multip, mesh);
    //Nextsim::VTK::write_dg("ResultsMEB_CGvel/Lambda", 0, Lambda, mesh);

    Nextsim::GlobalTimer.stop("time loop - i/o");

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
            std::cout << "\rAdvection step " << timestep << "\t " << std::setprecision(2)
                      << std::fixed << std::setw(10) << std::right << time << "s\t" << std::setw(8)
                      << std::right << timeInMinutes << "m\t" << std::setw(6) << std::right
                      << timeInHours << "h\t" << std::setw(6) << std::right << timeInDays << "d\t\t"
                      << std::flush;

        //! Initialize time-dependent forcing
        Nextsim::GlobalTimer.start("time loop - forcing");
        AtmForcingX.settime(time);
        AtmForcingY.settime(time);
        Nextsim::InterpolateCG(mesh, AX, AtmForcingX);
        Nextsim::InterpolateCG(mesh, AY, AtmForcingY);
        Nextsim::GlobalTimer.stop("time loop - forcing");

        Nextsim::GlobalTimer.start("time loop - advection");
        momentum.ProjectCGToDG(mesh, dgvx, vx);
        momentum.ProjectCGToDG(mesh, dgvy, vy);
        dgtransport.reinitvelocity();
        //dgtransport.step(dt_adv, A);
        dgtransport.step(dt_adv, H);
        dgtransport.step(dt_adv, D);

        //Nextsim::LimitMax(A, 1.0);
        //Nextsim::LimitMin(A, 0.0);
        Nextsim::LimitMin(H, 0.0);
        Nextsim::LimitMax(D, 1.0);
        Nextsim::LimitMin(D, 0.0);

        momentum.InterpolateDGToCG(mesh, cg_A, A);
        momentum.InterpolateDGToCG(mesh, cg_H, H);

        cg_A = cg_A.cwiseMin(1.0);
        cg_A = cg_A.cwiseMax(0.0);
        cg_H = cg_H.cwiseMax(1.e-4); //!< Limit H from below

        Nextsim::GlobalTimer.stop("time loop - advection");

        Nextsim::GlobalTimer.start("time loop - meb");

        //! Momentum subcycling
        for (size_t mom_step = 0; mom_step < mom_substeps; ++mom_step) {

            Nextsim::LimitMax(A, 1.0);
            Nextsim::LimitMin(A, 0.0);
            Nextsim::LimitMin(H, 0.0);
            Nextsim::LimitMax(D, 1.0);
            Nextsim::LimitMin(D, 0.0);

            Nextsim::GlobalTimer.start("time loop - meb - strain");
            //! Compute Strain Rate
            momentum.ProjectCG2VelocityToDG1Strain(mesh, E11, E12, E22, vx, vy);
            Nextsim::GlobalTimer.stop("time loop - meb - strain");

            Nextsim::GlobalTimer.start("time loop - meb - stress");
            //! Stress Update

            //Nextsim::MEB::StressUpdate(mesh, S11, S12, S22, E11, E12, E22,
            //    H, A, D, dt_momentum);

            //Nextsim::MEB::StressUpdateSandbox(mesh, S11, S12, S22, E11, E12, E22, H, A, D, DELTA,
            //    SHEAR, S1, S2, eta1, eta2, stressrelax, sigma_outside, tildeP, Pmax, td, d_crit,
            //    Regime, Multip, Lambda, dt_momentum);

            //Nextsim::MEB::StressUpdateMEB(mesh, S11, S12, S22, E11, E12, E22, H, A, D, DELTA,
            //    SHEAR, S1, S2, eta1, eta2, stressrelax, sigma_outside, tildeP, Pmax, td, d_crit,
            //    Regime, Multip, Lambda, dt_momentum);

            Nextsim::MEBSandbox::StressUpdateVP(mesh, S11, S12, S22, E11, E12, E22,
                H, A, RefScale::Pstar, RefScale::DeltaMin, dt_momentum);

            //Nextsim::MEB::StressUpdateMEB(mesh, S11, S12, S22, E11, E12, E22,
            //    H, A, D, dt_momentum);
            //Nextsim::MEB::ElasticUpdate(mesh, S11, S12, S22, E11, E12, E22,
            //    H, A, D, dt_momentum);

            Nextsim::GlobalTimer.stop("time loop - meb - stress");

            Nextsim::GlobalTimer.start("time loop - meb - update");
            //! Update
            Nextsim::GlobalTimer.start("time loop - meb - update1");

            //	    update by a loop.. explicit parts and h-dependent
            vx = (1.0
                / (RefScale::rho_ice * cg_H.array() / dt_momentum // implicit parts
                    + cg_A.array() * RefScale::F_ocean
                        * (OX.array() - vx.array()).abs()) // implicit parts
                * (RefScale::rho_ice * cg_H.array() / dt_momentum * vx.array() + //
                    cg_A.array()
                        * (RefScale::F_atm * AX.array().abs() * AX.array() + // atm forcing
                            RefScale::F_ocean * (OX - vx).array().abs()
                                * OX.array()) // ocean forcing
                    + RefScale::rho_ice * cg_H.array() * RefScaleCanada::fc
                        * (vx - OX).array() // cor + surface
                    ))
                     .matrix();
            vy = (1.0
                / (RefScale::rho_ice * cg_H.array() / dt_momentum // implicit parts
                    + cg_A.array() * RefScale::F_ocean
                        * (OY.array() - vy.array()).abs()) // implicit parts
                * (RefScale::rho_ice * cg_H.array() / dt_momentum * vy.array()
                    + cg_A.array()
                        * (RefScale::F_atm * AY.array().abs() * AY.array() + // atm forcing
                            RefScale::F_ocean * (OY - vy).array().abs()
                                * OY.array()) // ocean forcing
                    + RefScale::rho_ice * cg_H.array() * RefScaleCanada::fc
                        * (OY - vy).array() // cor + surface
                    ))
                     .matrix();
            Nextsim::GlobalTimer.stop("time loop - meb - update1");

            Nextsim::GlobalTimer.start("time loop - meb - update2");

            tmpx.zero();
            tmpy.zero();
            momentum.AddStressTensor(mesh, -1.0, tmpx, tmpy, S11, S12, S22);

            vx += (1.0
                / (RefScale::rho_ice * cg_H.array() / dt_momentum // implicit parts
                    + cg_A.array() * RefScale::F_ocean
                        * (OX.array() - vx.array()).abs()) // implicit parts
                * tmpx.array())
                      .matrix();
            vy += (1.0
                / (RefScale::rho_ice * cg_H.array() / dt_momentum // implicit parts
                    + cg_A.array() * RefScale::F_ocean
                        * (OY.array() - vy.array()).abs()) // implicit parts
                * tmpy.array())
                      .matrix();

            Nextsim::GlobalTimer.stop("time loop - meb - update2");
            Nextsim::GlobalTimer.stop("time loop - meb - update");

            Nextsim::GlobalTimer.start("time loop - meb - bound.");
            momentum.DirichletCompressionBottom(mesh, vx, 0.0);
            momentum.DirichletCompressionBottom(mesh, vy, 0.0);
            //momentum.DirichletCompressionFixCorner(mesh, vx); // Does not work
            //! inflow from top of concentrated ice
            momentum.DirichletCompressionTop(mesh, cg_H, 1.0);
            momentum.DirichletCompressionTop(mesh, H, 1.0);
            //momentum.DirichletCompressionTop(mesh, cg_A, 1.0);
            //momentum.DirichletCompressionTop(mesh, A, 1.0);
            momentum.DirichletCompressionTop(mesh, D, 0.0);

            Nextsim::GlobalTimer.stop("time loop - meb - bound.");
        }
        Nextsim::GlobalTimer.stop("time loop - meb");

        //OutputNextsim::Tools::Delta(mesh, E11, E12, E22, RefScaleCanada::DeltaMin, DELTA);
        Nextsim::VTK::write_dg("ResultsCompression/Delta", 0, DELTA, mesh);

        if (WRITE_VTK)
            if ((timestep % NT_vtk == 0)) {
                std::cout << "VTK output at day " << time / 24. / 60. / 60. << std::endl;

                size_t printstep = timestep / NT_vtk + 1.e-4;

                char s[80];
                sprintf(s, "ResultsCompression/ellipse_%03ld.txt", printstep);
                std::ofstream ELLOUT(s);
                for (size_t i = 0; i < mesh.n; ++i)
                    ELLOUT << Sinv1(i, 0) << "\t" << Sinv2(i, 0) << std::endl;
                ELLOUT.close();

                Nextsim::GlobalTimer.start("time loop - i/o");
                // Model Variables
                Nextsim::VTK::write_cg("ResultsCompression/vx", printstep, vx, mesh);
                Nextsim::VTK::write_cg("ResultsCompression/vy", printstep, vy, mesh);

                Nextsim::VTK::write_dg("ResultsCompression/A", printstep, A, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/H", printstep, H, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/D", printstep, D, mesh);

                Nextsim::VTK::write_dg("ResultsCompression/S11", printstep, S11, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/S12", printstep, S12, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/S22", printstep, S22, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/E11", printstep, E11, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/E12", printstep, E12, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/E22", printstep, E22, mesh);

                //save checkpoint
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/vx", printstep, vx);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/vy", printstep, vy);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/A", printstep, A);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/H", printstep, H);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/D", printstep, D);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/S11", printstep, S11);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/S12", printstep, S12);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/S22", printstep, S22);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/E11", printstep, E11);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/E12", printstep, E12);
                Nextsim::CheckPoints::saveData("ResultsCompression/Checkpoints/E22", printstep, E22);

                // Output Variables
                //Nextsim::Tools::Delta(mesh, E11, E12, E22, RefScaleCanada::DeltaMin, DELTA);
                //Nextsim::VTK::write_dg("ResultsCompression/Delta", printstep, DELTA, mesh);

                Nextsim::Tools::TensorInvariants(mesh, E11, E12, E22, Einv1, Einv2);
                Nextsim::VTK::write_dg("ResultsCompression/Einv1", printstep, Einv1, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/Einv2", printstep, Einv2, mesh);

                Nextsim::Tools::TensorInvariants(mesh, S11, S12, S22, Sinv1, Sinv2);
                Nextsim::VTK::write_dg("ResultsCompression/sigma_n", printstep, Sinv1, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/tau", printstep, Sinv2, mesh);

                //Nextsim::Tools::ElastoParams(
                //    mesh, E11, E12, E22, H, A, RefScaleCanada::DeltaMin, RefScaleCanada::Pstar, MU1, MU2);
                //Nextsim::VTK::write_dg("ResultsCompression/mu1", printstep, MU1, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/mu2", printstep, MU2, mesh);

                //Nextsim::VTK::write_dg(
                //    "ResultsCompression/stressOutside", printstep, sigma_outside, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/tildeP", printstep, tildeP, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/Pmax", printstep, Pmax, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/td", printstep, td, mesh);
                Nextsim::VTK::write_dg("ResultsCompression/d_crit", printstep, d_crit, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/eta1", printstep, eta1, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/eta2", printstep, eta2, mesh);
                //Nextsim::VTK::write_dg(
                //    "ResultsCompression/stressrelax", printstep, stressrelax, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/Regime", printstep, Regime, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/Multip", printstep, Multip, mesh);
                //Nextsim::VTK::write_dg("ResultsCompression/Lambda", printstep, Lambda, mesh);

                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
