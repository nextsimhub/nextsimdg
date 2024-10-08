/*!
 * @file AdvectionPeriodicBC_test.cpp
 * @date 07 Oct 2024
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "DGTransport.hpp"
#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "Tools.hpp"
#include "dgLimit.hpp"
#include "dgVisu.hpp"

#include <filesystem> // only for automatic creation of output directory
#include <iostream>

using namespace doctest;

/*!
 * Advection test case on a ring:
 * - center (0,0)
 * - inner radius 100 000 m
 * - outer radius 250 000 m
 *
 * Dirichlet conditions on inner and outer ring, periodic along ring
 *
 * Mesh is split into Nx elements in phi-direction and Ny elements in radial
 *
 *
 */

namespace ProblemConfig {
double R0 = 100000.0;
double R1 = 250000.0;

const size_t Nx0 = 32;
const size_t Ny0 = 4;
size_t Nx = Nx0;
size_t Ny = Ny0;

const size_t NT0 = 50;
size_t NT = NT0;
}

bool WRITE_VTK = true; //!< set to true for vtk output

double TOL = 1.e-7; //!< tolerance for checking test results

/*!
 *  Description of the test case
 *
 * circular transport in a domain of size [0,Lx] x [0,Ly]
 * Lx = 409600
 * Ly = 512000
 * (not using a square to test that non-square elements are handled correctly)
 *
 * Mesh of size Nx x Ny
 * Nx = 26, 52, 104
 * Ny = 32, 64, 128
 * (elements are not square to test this is working correctly)
 *
 * Time mesh
 * NT = 1000
 */

//! Packman-initial at 256000, 256000 with radius 128000
//! plus a smooth initial at 768000, 256000 with smaller radius
class Packman : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
        //      return 1;

        double r = (pow((x + 175000.0) / 50000, 2.0) + pow((y) / 50000.0, 2.0));
        if (r < 1.0)
            return exp(1.0) * exp(-1.0 / (1.0 - r));

        r = (pow((x) / 50000, 2.0) + pow((y - 175000.0) / 50000.0, 2.0));
        if (r < 1.0) {
            if (y > 175000)
                return 1.0;
            if (fabs(x) > 15000.0)
                return 1.0;
            return 0.0;
        }

        r = sqrt(pow((x - 175000.0) / 50000, 2.0) + pow((y) / 50000.0, 2.0));
        if (r < 1.0)
            return 1.0 - r;

        r = sqrt(pow((x) / 50000, 2.0) + pow((y + 175000.) / 50000.0, 2.0));
        if (r < 1.0)
            return 1.0;

        else
            return 0.0;
    }
};

// Velocity
class InitialVX : public Nextsim::Interpolations::Function { // (0.5,0.2) m/s

public:
    double operator()(double x, double y) const { return y * 2.0 * M_PI / ProblemConfig::R1; }
};
class InitialVY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const { return -x * 2.0 * M_PI / ProblemConfig::R1; }
};

//////////////////////////////////////////////////

template <int DG> class Test {

    const Nextsim::ParametricMesh& smesh; //!< Stores a reference to the spacial mesh.

    size_t N; //!< size of mesh N x N

    double dt; //!< time step size

    //! Velocity vectors and density
    Nextsim::DGVector<DG> phi;

    //! Transport main class
    Nextsim::DGTransport<DG> dgtransport;

    //! Velocity Field
    InitialVX VX;
    InitialVY VY;

    size_t writestep; //! write out n step in total (for debugging only)

public:
    Test(const Nextsim::ParametricMesh& mesh)
        : smesh(mesh)
        , dgtransport(smesh)
        , writestep(40)
    {
        //! Set time stepping scheme. 2nd order for dg0 and dg1, 3rd order dG2
        if (DG < 3)
            dgtransport.settimesteppingscheme("rk2");
        else
            dgtransport.settimesteppingscheme("rk3");
    }

    Test() { }

    void init()
    {
        dt = ProblemConfig::R1 / ProblemConfig::NT;
        //! Init Vectors
        phi.resize_by_mesh(smesh);
    }

    double run()
    {

        // init the test case, in particular resize vectors
        init();

        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, Packman());
        Nextsim::LimitMax(phi, 1.0);
        //	Nextsim::LimitMin(phi, 0.0);

        double initialmass = Nextsim::Tools::MeanValue(smesh, phi);
        // velocity field
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);

        std::string resultsdir = "test2_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);

        if (WRITE_VTK) {
            //! Compose name of output directory and create it
            std::filesystem::create_directory(resultsdir);

            Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", 0, phi, smesh);
        }
        std::cout << DG << "\t" << ProblemConfig::NT << "\t" << smesh.nx << "\t" << std::flush;

        // time loop
        for (size_t iter = 1; iter <= ProblemConfig::NT; ++iter) {

            dgtransport.reinitnormalvelocity();
            dgtransport.step(
                dt, phi); // performs one time step with the 2nd or 3rd Order Heun scheme
            Nextsim::LimitMax(phi, 1.0);
            Nextsim::LimitMin(phi, 0.0);

            if (WRITE_VTK)
                if (iter % (ProblemConfig::NT / writestep) == 0)
                    Nextsim::VTK::write_dg<DG>(
                        resultsdir + "/dg", iter / (ProblemConfig::NT / writestep), phi, smesh);
        }
        initialmass = (initialmass - Nextsim::Tools::MeanValue(smesh, phi)) / initialmass;

        std::cout << initialmass << "\t";
        // compute error
        return sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, Packman()))
            / ProblemConfig::R0;
    }
};

//////////////////////////////////////////////////

void create_rectanglemesh(Nextsim::ParametricMesh& smesh)
{
    smesh.reset(); // reset mesh and clear all variables

    // set number of elements and nodes
    smesh.nx = ProblemConfig::Nx;
    smesh.ny = ProblemConfig::Ny;
    smesh.nnodes = (ProblemConfig::Nx + 1) * (ProblemConfig::Ny + 1);
    smesh.nelements = ProblemConfig::Nx * ProblemConfig::Ny;

    // set vertices
    smesh.vertices.resize(smesh.nnodes, 2);
    size_t ii = 0;
    for (size_t iy = 0; iy <= ProblemConfig::Ny; ++iy)
        for (size_t ix = 0; ix <= ProblemConfig::Nx; ++ix, ++ii) {
            double r = ProblemConfig::R0
                + (ProblemConfig::R1 - ProblemConfig::R0) * iy / ProblemConfig::Ny;
            double p = -2.0 * M_PI * ix / ProblemConfig::Nx;
            smesh.vertices(ii, 0) = r * cos(p);
            smesh.vertices(ii, 1) = r * sin(p);
        }

    // landmask: set all to ice
    smesh.landmask.resize(smesh.nelements, 1);

    //// Boundary. Dirichlet on bottom / top, periodic left/right
    smesh.dirichlet[0].resize(smesh.nx); // bottom boundary
    for (size_t i = 0; i < smesh.nx; ++i)
        smesh.dirichlet[0][i] = i;

    smesh.dirichlet[2].resize(smesh.nx); // top boundary
    for (size_t i = 0; i < smesh.nx; ++i)
        smesh.dirichlet[2][i] = smesh.nx * (smesh.ny - 1) + i;

    smesh.periodic.resize(1);
    smesh.periodic[0].resize(smesh.ny);
    for (size_t i = 0; i < smesh.ny; ++i) {
        smesh.periodic[0][i][0] = 1; // periodic on left/right boundary
        smesh.periodic[0][i][1] = (i + 1) * smesh.nx
            - 1; // element on the left of the edge (so the most right one in the row)
        smesh.periodic[0][i][2]
            = i * smesh.nx; // element right of the edge (so the very left one) as...
        smesh.periodic[0][i][3]
            = i * (smesh.nx + 1) + i; // the edge is the most right edge in the row
    }
}

template <int DG> void run(const std::array<std::array<double, 6>, 3>& exact)
{
    Nextsim::ParametricMesh smesh(Nextsim::CARTESIAN);

#define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2 : -1))

    for (int it = 0; it < 2; ++it) {

        ProblemConfig::Nx = ProblemConfig::Nx0 * (1 << it);
        ProblemConfig::Ny = ProblemConfig::Ny0 * (1 << it);
        ProblemConfig::NT = ProblemConfig::NT0 * (DG2DEG(DG) + 1) * (DG2DEG(DG) + 1) * (1 << it);

        create_rectanglemesh(smesh);

        Test<DG> test(smesh);
        double error = test.run();
        std::cout << error << "\t" << exact[DG2DEG(DG)][it] << std::endl;
        CHECK(error == Approx(exact[DG2DEG(DG)][it]).epsilon(TOL));
    }
}

TEST_SUITE_BEGIN("Advection Periodic Boundary Conditions");
TEST_CASE("Advection Periodic Boundary Conditions")
{
    std::array<std::array<double, 6>, 3> exact = // Exact values taken 26/06/2024
        { std::array<double, 6>(
              { 1.0338503986019776e+00, 1.1451366598186583e+00, 1.0681593193338459e+00,
                  9.5252231195653603e-01, 8.1458581892610615e-01, 6.8950068528265651e-01 }),
            std::array<double, 6>(
                { 1.0949680727313791e+00, 8.3851748134278226e-01, 5.8501741891364523e-01,
                    4.0006092999256893e-01, 2.8493698219799068e-01, 2.0801424342840474e-01 }),
            std::array<double, 6>(
                { 7.1704413076190610e-01, 4.6135258391583606e-01, 3.2770311091970727e-01,
                    2.3424633076579465e-01, 1.7038068017215552e-01, 1.2486932724366147e-01 }) };

    std::cout << "DG\tNT\tNX\tmass\t\terror\t\texact" << std::endl;
    std::cout << std::setprecision(4) << std::scientific;
    run<1>(exact);
    run<3>(exact);
    run<6>(exact);
}
