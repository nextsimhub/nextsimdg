/*!
 * @file benchmark_box_mevp.cpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Tools.hpp"
#include "cgMomentum.hpp"
#include "cgVector.hpp"
#include "dgInitial.hpp"
#include "dgTransport.hpp"
#include "dgVisu.hpp"
#include "dgLimit.hpp"
#include "mevp.hpp"
#include "stopwatch.hpp"
#include "VPParameters.hpp"

#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

bool WRITE_VTK = true;
 
#define CG 2
#define DGadvection 1
#define DGstress 8

#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))

namespace Nextsim {
extern Timer GlobalTimer;
}

namespace ReferenceScale {
// Box-Test case from [Geosci. Model Dev., 8, 1747–1761, 2015]
constexpr double L = 1000000.0; //!< Size of domain
constexpr double vmax_ocean = 0.1; //!< Maximum velocity of ocean

constexpr double T = 30. * 24. * 60. * 60; //!< time

}

inline constexpr double SQR(double x) { return x * x; }
//! Description of the problem data, wind & ocean fields
struct OceanX {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (2.0 * y / ReferenceScale::L - 1);
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
        double X = M_PI * x / ReferenceScale::L;
        double Y = M_PI * y / ReferenceScale::L;
        constexpr double T = 4.0 * 24.0 * 60.0 * 60.0; //!< 4 days
        return 5.0 + (sin(2 * M_PI * time / T) - 3.0) * sin(2 * X) * sin(Y);
    }
};
struct AtmY {
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

struct InitialH {
public:
    double operator()(double x, double y) const
    {
        return 2.0; //!< constant ice height for box test
    }
};
struct InitialA {
public:
    double operator()(double x, double y) const { return x / ReferenceScale::L; }
};

int main()
{
    //!

    Nextsim::CGMomentum momentum;

    //! Define the spatial mesh
    Nextsim::Mesh mesh;
    constexpr size_t N = 32; //!< Number of mesh nodes
    mesh.BasicInit(N, N, ReferenceScale::L / N, ReferenceScale::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    //! define the time mesh
    constexpr double dt_adv = 120.0; //!< Time step of advection problem
    constexpr size_t NT = ReferenceScale::T / dt_adv + 1.e-4; //!< Number of Advections steps

    //! MEVP parameters
    constexpr double alpha = 300.0;
    constexpr double beta = 300.0;
    constexpr size_t NT_evp = 100;

    Nextsim::VPParameters VP;
    VP.C_atm = 2.25e-3; //!< See S. Danilov et al.: Finite element sea ice model. GMD 8, 2015
    VP.F_atm = VP.C_atm * VP.rho_atm;

    std::cout << "Time step size (advection) " << dt_adv << "\t" << NT << " time steps" << std::endl
              << "MEVP subcycling NTevp " << NT_evp << "\t alpha/beta " << alpha << " / " << beta
              << std::endl;

    //! VTK output
    constexpr double T_vtk = 1.0 * 24.0 * 60.0 * 60.0; // evey day
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

    // save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_cg("ResultsBox/vx", 0, vx, mesh);
    Nextsim::VTK::write_cg("ResultsBox/vy", 0, vy, mesh);
    Nextsim::VTK::write_dg("ResultsBox/A", 0, A, mesh);
    Nextsim::VTK::write_dg("ResultsBox/H", 0, H, mesh);
    Nextsim::VTK::write_dg("ResultsBox/Shear", 0, Nextsim::Tools::Shear(mesh, E11, E12, E22), mesh);
    Nextsim::VTK::write_dg("ResultsBox/Delta", 0, Nextsim::Tools::Delta(mesh, E11, E12, E22, VP.DeltaMin), mesh);
    Nextsim::GlobalTimer.stop("time loop - i/o");

    //! Transport
    Nextsim::CellVector<DGadvection> dgvx(mesh), dgvy(mesh);
    Nextsim::DGTransport<DGadvection, EDGEDOFS(DGadvection)> dgtransport(dgvx, dgvy);
    dgtransport.settimesteppingscheme("rk3");
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

        //! Initialize time-dependent forcing
        Nextsim::GlobalTimer.start("time loop - forcing");
        AtmForcingX.settime(time);
        AtmForcingY.settime(time);
        Nextsim::InterpolateCG(mesh, AX, AtmForcingX);
        Nextsim::InterpolateCG(mesh, AY, AtmForcingY);
        Nextsim::GlobalTimer.stop("time loop - forcing");

	//////////////////////////////////////////////////
        Nextsim::GlobalTimer.start("time loop - advection");
        momentum.ProjectCGToDG(mesh, dgvx, vx);
        momentum.ProjectCGToDG(mesh, dgvy, vy);
        dgtransport.reinitvelocity();
        dgtransport.step(dt_adv,A);
        dgtransport.step(dt_adv,H);

        A.col(0) = A.col(0).cwiseMin(1.0);
        A.col(0) = A.col(0).cwiseMax(0.0);
        H.col(0) = H.col(0).cwiseMax(0.0);

	//! Gauss-point limiting
        Nextsim::LimitMax(A, 1.0);
        Nextsim::LimitMin(A, 0.0);
        Nextsim::LimitMin(H, 0.0);
	
        momentum.InterpolateDGToCG(mesh, cg_A, A);
        momentum.InterpolateDGToCG(mesh, cg_H, H);
        cg_A = cg_A.cwiseMin(1.0);
        cg_A = cg_A.cwiseMax(0.0);
        cg_H = cg_H.cwiseMax(1.e-4); //!< Limit H from below

        Nextsim::GlobalTimer.stop("time loop - advection");	
	//////////////////////////////////////////////////

        Nextsim::GlobalTimer.start("time loop - mevp");
        //! Store last velocity for MEVP
        vx_mevp = vx;
        vy_mevp = vy;

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

            Nextsim::mEVP::StressUpdateHighOrder(VP, mesh, S11, S12, S22, E11, E12, E22, H, A, alpha, beta);

            Nextsim::GlobalTimer.stop("time loop - mevp - stress");

            Nextsim::GlobalTimer.start("time loop - mevp - update");
            //! Update
            Nextsim::GlobalTimer.start("time loop - mevp - update1");

            //	    update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
            for (int i = 0; i < vx.rows(); ++i) {
                vx(i) = (1.0
                    / (VP.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                        + cg_A(i) * VP.F_ocean
                            * fabs(OX(i) - vx(i))) // implicit parts
                    * (VP.rho_ice * cg_H(i) / dt_adv
                            * (beta * vx(i) + vx_mevp(i))
                        + // pseudo-timestepping
                        cg_A(i)
                            * (VP.F_atm * fabs(AX(i)) * AX(i) + // atm forcing
                                VP.F_ocean * fabs(OX(i) - vx(i))
                                    * OX(i)) // ocean forcing
                        + VP.rho_ice * cg_H(i) * VP.fc
                            * (vy(i) - OY(i)) // cor + surface
                        ));
                vy(i) = (1.0
                    / (VP.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                        + cg_A(i) * VP.F_ocean
                            * fabs(OY(i) - vy(i))) // implicit parts
                    * (VP.rho_ice * cg_H(i) / dt_adv
                            * (beta * vy(i) + vy_mevp(i))
                        + // pseudo-timestepping
                        cg_A(i)
                            * (VP.F_atm * fabs(AY(i)) * AY(i) + // atm forcing
                                VP.F_ocean * fabs(OY(i) - vy(i))
                                    * OY(i)) // ocean forcing
                        + VP.rho_ice * cg_H(i) * VP.fc
                            * (OX(i) - vx(i))));
            }
            Nextsim::GlobalTimer.stop("time loop - mevp - update1");

            Nextsim::GlobalTimer.start("time loop - mevp - update2");
            // Implicit etwas ineffizient
#pragma omp parallel for
            for (int i = 0; i < tmpx.rows(); ++i)
                tmpx(i) = tmpy(i) = 0;

            Nextsim::GlobalTimer.start("time loop - mevp - update2 -stress");
            momentum.AddStressTensor(mesh, -1.0, tmpx, tmpy, S11, S12, S22);
            Nextsim::GlobalTimer.stop("time loop - mevp - update2 -stress");

#pragma omp parallel for
            for (int i = 0; i < vx.rows(); ++i) {
                vx(i) += (1.0
                    / (VP.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                        + cg_A(i) * VP.F_ocean
                            * fabs(OX(i) - vx(i))) // implicit parts
                    * tmpx(i));

                vy(i) += (1.0
                    / (VP.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                        + cg_A(i) * VP.F_ocean
                            * fabs(OY(i) - vy(i))) // implicit parts
                    * tmpy(i));
            }
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

                int printstep = timestep / NT_vtk + 1.e-4;

                // char s[80];
                // sprintf(s, "ResultsBox/ellipse_%03d.txt", printstep);
                // std::ofstream ELLOUT(s);
                // for (size_t i = 0; i < mesh.n; ++i)
                //     ELLOUT << S1(i, 0) << "\t" << S2(i, 0) << std::endl;
                // ELLOUT.close();

                Nextsim::GlobalTimer.start("time loop - i/o");
                Nextsim::VTK::write_cg("ResultsBox/vx", printstep, vx, mesh);
                Nextsim::VTK::write_cg("ResultsBox/vy", printstep, vy, mesh);


                Nextsim::VTK::write_dg("ResultsBox/A", printstep, A, mesh);
                Nextsim::VTK::write_dg("ResultsBox/H", printstep, H, mesh);

		Nextsim::VTK::write_dg("ResultsBox/Shear", printstep, Nextsim::Tools::Shear(mesh, E11, E12, E22), mesh);
		Nextsim::VTK::write_dg("ResultsBox/Delta", printstep, Nextsim::Tools::Delta(mesh, E11, E12, E22, VP.DeltaMin), mesh);

                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
