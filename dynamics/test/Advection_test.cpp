/*!
 * @file Advection_test.cpp
 * @date 27 Aug 2024
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "DGTransport.hpp"
#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "Tools.hpp"
#include "dgLimiters.hpp"
#include "dgVisu.hpp"

#include <filesystem> // only for automatic creation of output directory
#include <iostream>
#include <vector>

using namespace doctest;

bool WRITE_VTK = false; //!< set to true for vtk output
int WRITE_EVERY = 5;

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
namespace ProblemConfig {
double Lx = 409600.0;
double Ly = 512000.0;
size_t Nx = 12;
const size_t Nx0 = 12;
size_t Ny = 13;
const size_t Ny0 = 13;
size_t NT = 100;
const size_t NT0 = 100;
}

//! Packman-initial at 256000, 256000 with radius 128000
//! plus a smooth initial at 768000, 256000 with smaller radius
class PackmanPlus : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
        //      return 1;
        double r = sqrt(pow((x - 240000.0) / 64000.0, 2.0) + pow((y - 120000.0) / 64000.0, 2.0));
        double r1 = sqrt(pow((x - 272000.0) / 96000.0, 2.0) + pow((y - 368000.0) / 96000.0, 2.0));
        if (r < 1.0)
            if ((16000.0 - fabs(x - 240000.0) < 0.0) || (y > 120000.0 - 16000.0))
                return 1.0 + exp(-5.0 * r1 * r1);
        return exp(-50.0 * r1 * r1);
    }
};
class SmoothBump : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
        double X = x / ProblemConfig::Lx;
        double Y = y / ProblemConfig::Lx;
        double r = (pow((X - 0.25), 2.0) + pow((Y - 0.5), 2.0)) / 0.025;
        if (r < 1)
            return exp(-1.0 / (1.0 - r));
        else
            return 0.0;
    }
};

// Velocity
class InitialVX : public Nextsim::Interpolations::Function { // (0.5,0.2) m/s

public:
    double operator()(double x, double y) const
    {
        return (y - 0.5 * ProblemConfig::Lx) * 2.0 * M_PI / ProblemConfig::Lx;
    }
};
class InitialVY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return (0.5 * ProblemConfig::Lx - x) * 2.0 * M_PI / ProblemConfig::Lx;
    }
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
        dt = ProblemConfig::Lx / ProblemConfig::NT;
        //! Init Vectors
        phi.resize_by_mesh(smesh);
    }

    double run()
    {

        // init the test case, in particular resize vectors
        init();

        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, SmoothBump());

        double initialaverage = Nextsim::Tools::MeanValue(smesh, phi);

        // velocity field
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);

        // Write out initial solution
        std::string resultsdir = "test1_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);
        if (WRITE_VTK) {
            //! Compose name of output directory and create it
            std::filesystem::create_directory(resultsdir);

            Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", 0, phi, smesh);
            Nextsim::VTK::write_dg<DG>(resultsdir + "/vx", 0, dgtransport.GetVx(), smesh);
            Nextsim::VTK::write_dg<DG>(resultsdir + "/vy", 0, dgtransport.GetVy(), smesh);
        }
        std::cout << DG << "\t" << ProblemConfig::NT << "\t" << smesh.nx << "\t" << std::flush;

        //! time loop
        for (size_t iter = 1; iter <= ProblemConfig::NT; ++iter) {

            dgtransport.reinitnormalvelocity();
            dgtransport.step(
                dt, phi); // performs one time step with the 2nd or 3rd Order Heun scheme
            if (WRITE_VTK)
                if (iter % (ProblemConfig::NT / writestep) == 0)
                    Nextsim::VTK::write_dg<DG>(
                        resultsdir + "/dg", iter / (ProblemConfig::NT / writestep), phi, smesh);
        }

        initialaverage = (initialaverage - Nextsim::Tools::MeanValue(smesh, phi))
            / initialaverage; // compute error in mass
        std::cout << initialaverage << "\t";

        // integrate the error
        return sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, SmoothBump()))
            / ProblemConfig::Lx;
    }
};

//////////////////////////////////////////////////

void create_rectanglemesh(Nextsim::ParametricMesh& smesh, const double Lx, const double Ly,
    const size_t Nx, const size_t Ny, double distort = 0.0)
{
    smesh.reset(); // reset mesh and clear all variables

    // set number of elements and nodes
    smesh.nx = Nx;
    smesh.ny = Ny;
    smesh.nnodes = (Nx + 1) * (Ny + 1);
    smesh.nelements = Nx * Ny;

    // set vertices
    smesh.vertices.resize(smesh.nnodes, 2);
    size_t ii = 0;
    for (size_t iy = 0; iy <= Ny; ++iy)
        for (size_t ix = 0; ix <= Nx; ++ix, ++ii) {
            smesh.vertices(ii, 0)
                = Lx * ix / Nx + Lx * distort * sin(M_PI * ix / Nx * 3.0) * sin(M_PI * iy / Ny);
            smesh.vertices(ii, 1) = Ly * iy / Ny
                + Ly * distort * sin(M_PI * iy / Ny * 2.0) * sin(M_PI * ix / Nx * 2.0);
        }

    //// Boundary
    // landmask: set all to ice
    smesh.landmask.resize(smesh.nelements, 1);

    smesh.dirichlet[0].resize(smesh.nx); // bottom boundary
    for (size_t i = 0; i < smesh.nx; ++i)
        smesh.dirichlet[0][i] = i;

    smesh.dirichlet[1].resize(smesh.ny); // right boundary
    for (size_t i = 0; i < smesh.ny; ++i)
        smesh.dirichlet[1][i] = i * Nx + Nx - 1;

    smesh.dirichlet[2].resize(smesh.nx); // top boundary
    for (size_t i = 0; i < smesh.nx; ++i)
        smesh.dirichlet[2][i] = Nx * (Ny - 1) + i;

    smesh.dirichlet[3].resize(smesh.ny); // left boundary
    for (size_t i = 0; i < smesh.ny; ++i)
        smesh.dirichlet[3][i] = i * Nx;

    // no other boundaries.
}

template <int DG> void run(double distort, const std::array<std::array<double, 6>, 3>& exact)
{
    Nextsim::ParametricMesh smesh(Nextsim::CARTESIAN); // 0 means no output

#define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2 : -1))
    for (int it = 0; it < 2; ++it) {

        // set number of elements in x- and y-direction
        ProblemConfig::Nx = ProblemConfig::Nx0 * (1 << it);
        ProblemConfig::Ny = ProblemConfig::Ny0 * (1 << it);
        // set number of time-steps such that the cfl condition is fulfilled
        ProblemConfig::NT = ProblemConfig::NT0 * (DG2DEG(DG) + 1) * (DG2DEG(DG) + 1) * (1 << it);

        // create the mesh
        create_rectanglemesh(smesh, ProblemConfig::Lx, ProblemConfig::Ly, ProblemConfig::Nx,
            ProblemConfig::Ny, distort);

        Test<DG> test(smesh);
        double error = test.run();
        std::cout << error << "\t" << exact[DG2DEG(DG)][it] << std::endl;
        CHECK(error == Approx(exact[DG2DEG(DG)][it]).epsilon(TOL));
    }
}

TEST_SUITE_BEGIN("Advection");
TEST_CASE("Advection")
{
    std::cout << std::setprecision(4) << std::scientific;

    // Exact values taken 25/06/2024 after some plausibility check of the results
    // and by checking the theoretical order of convergence to be expected
    std::array<std::array<double, 6>, 3> exact
        = { std::array<double, 6>(
                { 5.1256149074257538e-02, 4.8288256703303903e-02, 4.3105635248809886e-02,
                    3.5793986049422598e-02, 2.6937487824223016e-02, 1.8456775604583933e-02 }),
              std::array<double, 6>(
                  { 2.9450967798560313e-02, 1.3427281824939470e-02, 5.8574889800512000e-03,
                      2.2756813704797301e-03, 7.3564079504566610e-04, 1.9177169540867740e-04 }),
              std::array<double, 6>(
                  { 9.9340386651904228e-03, 4.0274889287136816e-03, 1.2400186092664943e-03,
                      2.8460130790583933e-04, 4.5459555769938981e-05, 4.6851425365137733e-06 }) };

    std::cout << std::endl << "DG\tNT\tNX\tmass loss\terror\t\texact" << std::endl;
    run<1>(0.0, exact);
    run<3>(0.0, exact);
    run<6>(0.0, exact);
}
TEST_CASE("Distorted Mesh")
{
    std::array<std::array<double, 6>, 3> exact
        = { std::array<double, 6>(
                { 5.0958748236875594e-02, 4.8461087243594887e-02, 4.3918572621094686e-02,
                    3.7038043078562323e-02, 2.8334464207979679e-02, 1.9666196904181983e-02 }),
              std::array<double, 6>(
                  { 3.2033423984965226e-02, 1.5639052766007147e-02, 6.7926997203524983e-03,
                      2.7369916393700151e-03, 9.2109964949992139e-04, 2.5128064874203957e-04 }),
              std::array<double, 6>(
                  { 1.1830071142946669e-02, 4.9207680503709564e-03, 1.5831512504162868e-03,
                      3.8869595113351363e-04, 6.7162895595772516e-05, 7.5853786764529606e-06 }) };

    std::cout << std::endl << "Distorted mesh" << std::endl;
    std::cout << "DG\tNT\tNX\tmass loss\terror\t\texact" << std::endl;
    run<1>(0.05, exact);
    run<3>(0.05, exact);
    run<6>(0.05, exact);
}
