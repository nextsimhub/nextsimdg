/*!
 * @file benchmark_mehlmann_mevp.cpp
 * @date 21 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#include "Tools.hpp"
#include "cgMomentum.hpp"
#include "cgVector.hpp"
#include "dgInitial.hpp"
#include "dgTransport.hpp"
#include "dgVisu.hpp"
#include "mMEB.hpp"
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

inline constexpr double SQR(double x) { return x * x; }

//! Description of the problem data, wind & ocean fields
struct OceanX {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (2.0 * y / ReferenceScale::L - 1.0);
    }
};
struct OceanY {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (1.0 - 2.0 * x / ReferenceScale::L);
    }
};

struct AtmX {
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
struct AtmY {
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
struct InitialH {
public:
    double operator()(double x, double y) const
    {
        return 0.3 + 0.005 * (sin(6.e-5 * x) + sin(3.e-5 * y));
    }
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
    constexpr size_t N = 128; //!< Number of mesh nodes
    mesh.BasicInit(N, N, ReferenceScale::L / N, ReferenceScale::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    //! define the time mesh
    constexpr double dt_adv = 120.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    //constexpr double alpha = 48000.0;
    constexpr double alpha = 800.0;
    constexpr double beta = 800.0;
    constexpr size_t NT_evp = 1000;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta
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
    Nextsim::CGVector<CG> vx_mevp(mesh), vy_mevp(mesh); //!< temp. Velocity used for MEVP
    Nextsim::CGVector<CG> vx_p(mesh), vy_p(mesh);
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

    Nextsim::CellVector<DGadvection> D(mesh); //!< ice damage. ?? Really dG(0) ??
    Nextsim::L2ProjectInitial(mesh, D, InitialD());

    Nextsim::CellVector<DGstress> S11_mmeb(mesh), S12_mmeb(mesh), S22_mmeb(mesh);
    Nextsim::CellVector<DGstress> S11_p(mesh), S12_p(mesh), S22_p(mesh);
    S11.zero();
    S12.zero();
    S22.zero();
    S11_mmeb = S11;
    S12_mmeb = S12;
    S22_mmeb = S22;

    Nextsim::CellVector<1> DELTA(mesh); //!< Storing DELTA
    Nextsim::CellVector<1> SHEAR(mesh); //!< Storing DELTA
    Nextsim::CellVector<1> S1(mesh), S2(mesh); //!< Stress invariants
    Nextsim::CellVector<1> MU1(mesh), MU2(mesh); //!< Stress invariants

    // save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg("Results_mMEB/vx", 0, vx, mesh);
    Nextsim::VTK::write_cg("Results_mMEB/vy", 0, vy, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/A", 0, A, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/H", 0, H, mesh);

    Nextsim::Tools::Delta(mesh, E11, E12, E22, ReferenceScale::DeltaMin, DELTA);
    Nextsim::VTK::write_dg("Results_mMEB/Delta", 0, DELTA, mesh);
    Nextsim::Tools::Shear(mesh, E11, E12, E22, ReferenceScale::DeltaMin, SHEAR);
    Nextsim::VTK::write_dg("Results_mMEB/Shear", 0, SHEAR, mesh);
    Nextsim::Tools::ElastoParams(
        mesh, E11, E12, E22, H, A, ReferenceScale::DeltaMin, ReferenceScale::Pstar, MU1, MU2);
    Nextsim::VTK::write_dg("Results_mMEB/mu1", 0, MU1, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/mu2", 0, MU2, mesh);

    Nextsim::VTK::write_dg("Results_mMEB/S11", 0, S11, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/S12", 0, S12, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/S22", 0, S22, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/E11", 0, E11, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/E12", 0, E12, mesh);
    Nextsim::VTK::write_dg("Results_mMEB/E22", 0, E22, mesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    //! Transport
    Nextsim::CellVector<DGadvection> dgvx(mesh), dgvy(mesh);
    Nextsim::DGTransport<DGadvection, DGadvection> dgtransport(dgvx, dgvy);
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
            std::cout << "\rAdvection step " << timestep << "\t " << std::setprecision(2)
                      << std::fixed << std::setw(10) << std::right << time << "s\t" << std::setw(8)
                      << std::right << timeInMinutes << "m\t" << std::setw(6) << std::right
                      << timeInHours << "h\t" << std::setw(6) << std::right << timeInDays << "d\t\t"
                      << std::flush;

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
        dgtransport.step(dt_adv, D);

        //! Very simple limiting (just constants)
        A.col(0) = A.col(0).cwiseMin(1.0);
        A.col(0) = A.col(0).cwiseMax(0.0);
        H.col(0) = H.col(0).cwiseMax(0.0);
        D.col(0) = D.col(0).cwiseMin(1.0);
        D.col(0) = D.col(0).cwiseMax(0.0);

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

        //Nextsim::MEB::mMEBStressRelaxation(mesh, S11, S12, S22, D, A, dt_adv);
        S11_mmeb = S11;
        S12_mmeb = S12;
        S22_mmeb = S22;

        //! MEVP subcycling
        Nextsim::mMEB::mMEBStressRelaxation(mesh, S11, S12, S22, D, A, dt_adv);

        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {

            //Check now close are v^p and v^{p-1} and S^p S^{p-1}
            if ((mevpstep + 1) % NT_evp == 0) {
                std::cout << "Advection step " << timestep << " mEV iteration " << mevpstep + 1 << std::endl;
                std::cout << "Norm vx = " << std::setprecision(6) << (vx - vx_p).norm() << std::endl;
                std::cout << "Norm vy = " << (vy - vy_p).norm() << std::endl;
                std::cout << "Norm S11 = " << (S11 - S11_p).norm() << std::endl;
                std::cout << "Norm S12 = " << (S12 - S12_p).norm() << std::endl;
                std::cout << "Norm S22 = " << (S22 - S22_p).norm() << std::endl;
            }
            vx_p = vx;
            vy_p = vy;
            S11_p = S11;
            S12_p = S12;
            S22_p = S22;

            Nextsim::GlobalTimer.start("time loop - mevp - strain");
            //! Compute Strain Rate
            momentum.ProjectCG2VelocityToDG1Strain(mesh, E11, E12, E22, vx, vy);
            Nextsim::GlobalTimer.stop("time loop - mevp - strain");

            Nextsim::GlobalTimer.start("time loop - mevp - update");
            //! Update
            Nextsim::GlobalTimer.start("time loop - mevp - update1");

            //	    update by a loop.. implicit parts and h-dependent

            //
            vx = (1.0
                / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A.array() * ReferenceScale::F_ocean
                        * (OX.array() - vx.array()).abs()) // implicit parts
                * (ReferenceScale::rho_ice * cg_H.array() / dt_adv
                        * (beta * vx.array() + vx_mevp.array())
                    + // pseudo-timestepping
                    cg_A.array()
                        * (ReferenceScale::F_atm * AX.array().abs() * AX.array() + // atm forcing
                            ReferenceScale::F_ocean * (OX - vx).array().abs()
                                * OX.array()) // ocean forcing
                    + ReferenceScale::rho_ice * cg_H.array() * ReferenceScale::fc
                        * (vy - OY).array() // cor + surface
                    ))
                     .matrix();
            vy = (1.0
                / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A.array() * ReferenceScale::F_ocean
                        * (OY.array() - vy.array()).abs()) // implicit parts
                * (ReferenceScale::rho_ice * cg_H.array() / dt_adv
                        * (beta * vy.array() + vy_mevp.array())
                    + // pseudo-timestepping
                    cg_A.array()
                        * (ReferenceScale::F_atm * AY.array().abs() * AY.array() + // atm forcing
                            ReferenceScale::F_ocean * (OY - vy).array().abs()
                                * OY.array()) // ocean forcing
                    + ReferenceScale::rho_ice * cg_H.array() * ReferenceScale::fc
                        * (OX - vx).array() // cor + surface
                    ))
                     .matrix();
            Nextsim::GlobalTimer.stop("time loop - mevp - update1");

            Nextsim::GlobalTimer.start("time loop - mevp - update2");
            // Implicit etwas ineffizient
            tmpx.zero();
            tmpy.zero();
            momentum.AddStressTensor(mesh, -1.0, tmpx, tmpy, S11, S12, S22);
            vx += (1.0
                / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A.array() * ReferenceScale::F_ocean
                        * (OX.array() - vx.array()).abs()) // implicit parts
                * tmpx.array())
                      .matrix();
            vy += (1.0
                / (ReferenceScale::rho_ice * cg_H.array() / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A.array() * ReferenceScale::F_ocean
                        * (OY.array() - vy.array()).abs()) // implicit parts
                * tmpy.array())
                      .matrix();

            Nextsim::GlobalTimer.stop("time loop - mevp - update2");
            Nextsim::GlobalTimer.stop("time loop - mevp - update");

            Nextsim::GlobalTimer.start("time loop - mevp - stress");

            Nextsim::mMEB::mMEBStressUpdate(mesh, S11, S12, S22, E11, E12, E22, H, A, D,
                dt_adv, alpha, S11_mmeb, S12_mmeb, S22_mmeb);
            Nextsim::mMEB::mMEBDamageUpdate(mesh, S11, S12, S22, D, A, dt_adv);

            //Nextsim::MEB::mMEBStressRelaxation(mesh, S11, S12, S22, D, A, dt_adv);

            //Nextsim::mEVP::StressUpdate(mesh, S11, S12, S22, E11, E12, E22, H, A,
            //    ReferenceScale::Pstar, ReferenceScale::DeltaMin, alpha, beta);

            Nextsim::GlobalTimer.stop("time loop - mevp - stress");

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

                int printstep = timestep / NT_vtk + 1.e-4;

                char s[80];
                sprintf(s, "Results_mMEB/ellipse_%03d.txt", printstep);
                std::ofstream ELLOUT(s);
                for (size_t i = 0; i < mesh.n; ++i)
                    ELLOUT << S1(i, 0) << "\t" << S2(i, 0) << std::endl;
                ELLOUT.close();

                Nextsim::GlobalTimer.start("time loop - i/o");
                Nextsim::VTK::write_cg("Results_mMEB/vx", printstep, vx, mesh);
                Nextsim::VTK::write_cg("Results_mMEB/vy", printstep, vy, mesh);
                Nextsim::VTK::write_dg("Results_mMEB/A", printstep, A, mesh);
                Nextsim::VTK::write_dg("Results_mMEB/H", printstep, H, mesh);
                Nextsim::VTK::write_dg("Results_mMEB/D", printstep, D, mesh);

                Nextsim::Tools::Delta(mesh, E11, E12, E22, ReferenceScale::DeltaMin, DELTA);
                Nextsim::VTK::write_dg("Results_mMEB/Delta", printstep, DELTA, mesh);
                Nextsim::Tools::Shear(mesh, E11, E12, E22, ReferenceScale::DeltaMin, SHEAR);
                Nextsim::VTK::write_dg("Results_mMEB/Shear", printstep, SHEAR, mesh);

                // Nextsim::Tools::ElastoParams(mesh, E11, E12, E22, H, A,
                //     ReferenceScale::DeltaMin, ReferenceScale::Pstar, MU1, MU2);
                // Nextsim::VTK::write_dg("Results_mMEB/mu1", printstep, MU1, mesh);
                // Nextsim::VTK::write_dg("Results_mMEB/mu2", printstep, MU2, mesh);

                Nextsim::VTK::write_dg("Results_mMEB/S11", printstep, S11, mesh);
                Nextsim::VTK::write_dg("Results_mMEB/S12", printstep, S12, mesh);
                Nextsim::VTK::write_dg("Results_mMEB/S22", printstep, S22, mesh);
                // Nextsim::VTK::write_dg("Results_mMEB/E11", printstep, E11, mesh);
                // Nextsim::VTK::write_dg("Results_mMEB/E12", printstep, E12, mesh);
                // Nextsim::VTK::write_dg("Results_mMEB/E22", printstep, E22, mesh);

                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}