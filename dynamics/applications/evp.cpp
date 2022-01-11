#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "benchmark_data.hpp"
#include "cginitial.hpp"
#include "cgmomentum.hpp"
#include "cgvector.hpp"
#include "dgvisu.hpp"
#include "dynamics.hpp"
#include "stopwatch.hpp"

bool WRITE_VTK = true;

namespace Nextsim {
extern Timer GlobalTimer;
}

int main()
{
    Nextsim::Dynamics dynamics;

    constexpr size_t N = 100; //!< Number of mesh nodes
    dynamics.GetMesh().BasicInit(N, N, ReferenceScale::L / N, ReferenceScale::L / N);
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Spatial mesh with mesh " << N << " x " << N << " elements." << std::endl;

    constexpr double T = 2.0 * 24 * 60 * 60; //!< Time horizon (in sec)
    constexpr double k_adv = 90.0; //!< Time step of advection problem

    // Compute Parabolic CFL
    // get the effective viscosity, e.g. d_t v = 2 eta / rho_ice div(sigma)
    constexpr double eta = 1.e12;

    constexpr double effective_viscosity = 2.0 * eta / ReferenceScale::rho_ice;
    constexpr double gamma = 5; //!< Penalty parameter for internal continuity
    constexpr double gammaboundary = gamma; //!< Penalty parameter for boundary data
    constexpr double cfl = 0.05 / gamma / effective_viscosity;
    //! compute time step. Make the advection step a multiple of it.
    constexpr double k = k_adv / static_cast<int>(0.5 + k_adv / (SQR(ReferenceScale::L / N) * cfl));
    constexpr size_t NT = T / k + 1.e-4;
    constexpr size_t NTadv = k_adv / k + 1.e-4;

    std::cout << "Time step size " << k << "\t" << NT << " time steps" << std::endl
              << "Advection step " << k_adv << "\t every " << NTadv << " momentum steps"
              << std::endl;

    dynamics.GetTimeMesh().BasicInit(NT, k_adv, k);

    //! VTK output every hour
    constexpr double Tvtk = 1.0 * 60.0 * 60.0;
    constexpr size_t NTvtk = Tvtk / k;

    //! Initialize the Dynamical Core (vector sizes, etc.)
    dynamics.BasicInit();

    //! Initial data of the problem
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetH(), InitialH());
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetA(), InitialA());
    AtmX AtmForcingX; //!< stupid names....
    AtmY AtmForcingY;
    AtmForcingX.settime(0.0);
    AtmForcingY.settime(0.0);
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmX(), AtmForcingX);
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmY(), AtmForcingY);
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetOceanX(), OceanX());
    Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetOceanY(), OceanY());

    //! Initialize the velocity
    dynamics.GetVX().zero();
    dynamics.GetVY().zero();

    // save initial condition
    Nextsim::GlobalTimer.start("time loop - i/o");
    Nextsim::VTK::write_dg<2>("Results/vx", 0, dynamics.GetVX(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<2>("Results/vy", 0, dynamics.GetVY(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<2>("Results/A", 0, dynamics.GetA(), dynamics.GetMesh());
    Nextsim::VTK::write_dg<2>("Results/H", 0, dynamics.GetH(), dynamics.GetMesh());
    Nextsim::GlobalTimer.stop("time loop - i/o");

    Nextsim::GlobalTimer.start("time loop");
    size_t advectionstep = 0;

    for (size_t timestep = 1; timestep <= dynamics.GetTimeMesh().N; ++timestep) {
        Nextsim::GlobalTimer.start("time loop - reinit");
        double time = dynamics.GetTimeMesh().dt_momentum * timestep;
        if (timestep % 1000000 == 0)
            std::cout << "--- Time step " << timestep << "\t advection step " << advectionstep
                      << "-> day " << time / (24.0 * 60.0 * 60.0) << std::endl;

        // advection step?
        if (timestep % NTadv == 0) {
            ++advectionstep;

            //! Initial (atm) Forcing (ocean is stationary)
            AtmForcingX.settime(time);
            AtmForcingY.settime(time);
            Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmX(), AtmForcingX);
            Nextsim::L2ProjectInitial(dynamics.GetMesh(), dynamics.GetAtmY(), AtmForcingY);
            Nextsim::GlobalTimer.stop("time loop - reinit");

            // transport step
            dynamics.advectionStep();
        }

        ////// momentum step

        // set the rhs forcing
        dynamics.GetTMPX().zero();
        dynamics.GetTMPY().zero();

        // forcing
#pragma omp parallel for
        for (size_t i = 0; i < dynamics.GetMesh().n; ++i) {
            assert(dynamics.GetH()(i, 0) > 1.e-6);

            double scale = dynamics.GetA()(i, 0) / dynamics.GetH()(i, 0);

            dynamics.GetTMPX()(i, 0) = scale * ReferenceScale::C_atm * ReferenceScale::rho_atm
                / ReferenceScale::rho_ice * dynamics.GetAtmX()(i, 0)
                * abs(dynamics.GetAtmX()(i, 0));
            dynamics.GetTMPY()(i, 0) = scale * ReferenceScale::C_atm * ReferenceScale::rho_atm
                / ReferenceScale::rho_ice * dynamics.GetAtmY()(i, 0)
                * abs(dynamics.GetAtmY()(i, 0));

            dynamics.GetTMPX()(i, 0) += scale * ReferenceScale::C_ocean * ReferenceScale::rho_ocean
                / ReferenceScale::rho_ice * dynamics.GetOceanX()(i, 0)
                * abs(dynamics.GetOceanX()(i, 0) - dynamics.GetVX()(i, 0));
            dynamics.GetTMPY()(i, 0) += scale * ReferenceScale::C_ocean * ReferenceScale::rho_ocean
                / ReferenceScale::rho_ice * dynamics.GetOceanY()(i, 0)
                * abs(dynamics.GetOceanY()(i, 0) - dynamics.GetVY()(i, 0));

            for (size_t c = 0; c < 6; ++c) {
                dynamics.GetTMPX()(i, c) -= scale * ReferenceScale::C_ocean
                    * ReferenceScale::rho_ocean / ReferenceScale::rho_ice * dynamics.GetVX()(i, c)
                    * abs(dynamics.GetOceanX()(i, 0) - dynamics.GetVX()(i, 0));
                dynamics.GetTMPY()(i, c) -= scale * ReferenceScale::C_ocean
                    * ReferenceScale::rho_ocean / ReferenceScale::rho_ice * dynamics.GetVY()(i, c)
                    * abs(dynamics.GetOceanY()(i, 0) - dynamics.GetVY()(i, 0));
            }
        }

        //! compute stress
        dynamics.computeStrainRateTensor(); //!< E = 1/2 (nabla v + nabla v^T)

        //! compute viscosity
        Nextsim::CellVector<0> zeta(dynamics.GetMesh());
#pragma omp parallel for
        for (size_t i = 0; i < dynamics.GetMesh().n; ++i) {
            double delta = sqrt(0.5
                    * (0.5 * SQR(dynamics.GetE11()(i, 0) - dynamics.GetE22()(i, 0))
                        + 2.0 * SQR(dynamics.GetE12()(i, 0)))
                + SQR(dynamics.GetE11()(i, 0) + dynamics.GetE22()(i, 0)) + SQR(2.e-9));
            zeta(i, 0) = 0.5 / delta * ReferenceScale::Pstar * dynamics.GetH()(i, 0)
                * exp(-20.0 * (1.0 - dynamics.GetA()(i, 0)));
        }

        double stressscale = 1.0 / ReferenceScale::rho_ice;

        // compute stress

#pragma omp parallel for
        for (size_t i = 0; i < dynamics.GetMesh().n; ++i) {
            double Pi = 0.5 * ReferenceScale::Pstar * dynamics.GetH()(i, 0) * exp(-20.0 * (1.0 - dynamics.GetA()(i, 0)));
            for (size_t c = 0; c < 3; ++c) {
                dynamics.GetS11()(i, c) = zeta(i, 0) * (0.5 * dynamics.GetE11()(i, c) + 0.75 * (dynamics.GetE11()(i, c) + dynamics.GetE22()(i, c))) - Pi;
                dynamics.GetS12()(i, c) = zeta(i, 0) * (0.5 * dynamics.GetE12()(i, c));
                dynamics.GetS22()(i, c) = zeta(i, 0) * (0.5 * dynamics.GetE22()(i, c) + 0.75 * (dynamics.GetE11()(i, c) + dynamics.GetE22()(i, c))) - Pi;
            }
        }

        dynamics.addStressTensor(-1.0 * stressscale); //!< tmp += div(S)
        dynamics.velocityContinuity(eta * gamma * stressscale); //!< tmp += < [v], [phi] >
        dynamics.velocityDirichletBoundary(eta * gammaboundary * stressscale); //!< tmp += < v, [phi] >_G

        //! Velocity update by explicit Euler
        dynamics.GetVX() += dynamics.GetTimeMesh().dt_momentum * dynamics.GetTMPX();
        dynamics.GetVY() += dynamics.GetTimeMesh().dt_momentum * dynamics.GetTMPY();

        //! Output
        if (WRITE_VTK)
            if ((timestep % NTvtk == 0)) {

                size_t printstep = timestep / NTvtk;
                Nextsim::GlobalTimer.start("time loop - i/o");

                Nextsim::VTK::write_dg<2>(
                    "Results/vx", printstep, dynamics.GetVX(), dynamics.GetMesh());
                Nextsim::VTK::write_dg<2>(
                    "Results/vy", printstep, dynamics.GetVY(), dynamics.GetMesh());
                // Nextsim::VTK::write_dg<1>("Results/S11", printstep, dynamics.GetS11(),
                // dynamics.GetMesh()); Nextsim::VTK::write_dg<1>("Results/S12", printstep,
                // dynamics.GetS12(), dynamics.GetMesh()); Nextsim::VTK::write_dg<1>("Results/S22",
                // printstep, dynamics.GetS22(), dynamics.GetMesh());
                Nextsim::VTK::write_dg<2>(
                    "Results/A", printstep, dynamics.GetA(), dynamics.GetMesh());
                Nextsim::VTK::write_dg<2>(
                    "Results/H", printstep, dynamics.GetH(), dynamics.GetMesh());
                Nextsim::VTK::write_dg<0>("Results/zeta", printstep, zeta, dynamics.GetMesh());

                // Nextsim::VTK::write_dg<0>("Results/D",printstep,dynamics.GetD(),
                // dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/ox",printstep,dynamics.GetOceanX(),
                // dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/oy",printstep,dynamics.GetOceanY(),
                // dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/ax",printstep,dynamics.GetAtmX(),
                // dynamics.GetMesh());
                // Nextsim::VTK::write_dg<0>("Results/ay",printstep,dynamics.GetAtmY(),
                // dynamics.GetMesh());

                Nextsim::GlobalTimer.stop("time loop - i/o");
            }
    }
    Nextsim::GlobalTimer.stop("time loop");

    std::cout << std::endl;
    Nextsim::GlobalTimer.print();
}
