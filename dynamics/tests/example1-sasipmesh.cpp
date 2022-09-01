/*!
 * @file eexample1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "Mesh.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "ParametricTransport.hpp"
#include "dgLimiters.hpp"
#include "dgVisu.hpp"

#include "stopwatch.hpp"
#include "testtools.hpp"
#include <cassert>
#include <chrono>
#include <filesystem> // only for automatic creation of output directory
#include <iostream>
#include <vector>

namespace Nextsim {
extern Timer GlobalTimer;
}

bool WRITE_VTK = true; //!< set to true for vtk output
int WRITE_EVERY = 5;

double TOL = 1.e-10; //!< tolerance for checking test results

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
 */
namespace ProblemConfig {
double Lx = 409600.0;
double Ly = 512000.0;
size_t Nx = 26;
size_t Ny = 32;
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
        return exp(-25.0 * r1 * r1);
    }
};
class SmoothBump : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
        double r = (sqrt(pow((x - 0.25 * ProblemConfig::Lx), 2.0) + pow((y - 0.5 * ProblemConfig::Ly), 2.0))) / ProblemConfig::Lx / 0.2;
        if (r < 1)
            return exp(-1.0 / (1.0 - r));
        else
            return 0.0;
    }
};

// Velocity
class InitialVX : public Nextsim::Interpolations::Function { // (0.5,0.2) m/s

    double _time;

public:
    void settime(double t) { _time = t; }

    double operator()(double x, double y) const
    {
        return (y - 0.5 * ProblemConfig::Lx) * 2.0 * M_PI;
    }
};
class InitialVY : public Nextsim::Interpolations::Function {

    double _time;

public:
    void settime(double t) { _time = t; }

    double operator()(double x, double y) const
    {
        return (0.5 * ProblemConfig::Lx - x) * 2.0 * M_PI;
    }
};

//////////////////////////////////////////////////

template <int DG>
class Test {

    const Nextsim::ParametricMesh& smesh; //!< Stores a reference to the spacial mesh.

    size_t N; //!< size of mesh N x N

    size_t NT; //!< number of time steps
    double dt; //!< time step size

    Nextsim::Mesh mesh; //!< space mesh

    //! Velocity vectors and density
    Nextsim::CellVector<DG> phi;

    //! Transport main class
    Nextsim::ParametricTransport<DG> dgtransport;

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
        //! Init Time Mesh
        double cfl = 0.01;

        double hmin = smesh.hmin();

        dt = cfl * hmin / (0.5 * ProblemConfig::Lx);
        double tmax = 1.0;
        NT = (static_cast<int>((tmax / dt + 1) / 100 + 1) * 100); // No time steps dividable by 100
        dt = tmax / NT;
        //! Init Vectors
        phi.resize_by_mesh(smesh);
    }

    void run()
    {

        //! Compose name of output directory and create it
        std::string resultsdir = "Example1_" + std::to_string(DG) + "_" + std::to_string(mesh.nx);
        std::filesystem::create_directory(resultsdir);

        Nextsim::GlobalTimer.reset();

        // init the test case, in particular resize vectors
        init();

        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, SmoothBump());

        if (WRITE_VTK) {
            Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", 0, phi, smesh);
        }

        //! time loop
        for (size_t iter = 1; iter <= NT; ++iter) {
            VX.settime(dt * iter);
            VY.settime(dt * iter);
            Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
            Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);
            dgtransport.reinitnormalvelocity();

            dgtransport.step(dt, phi); // performs one time step with the 2nd or 3rd Order Heun scheme
            if (WRITE_VTK)
                if (iter % (NT / writestep) == 0)
                    Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", iter / (NT / writestep), phi, smesh);
        }

        // integrate the error
        std::cout << Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, SmoothBump()) << std::endl;
    }
};

//////////////////////////////////////////////////

void create_rectanglemesh(const double Lx, const double Ly, const size_t Nx, const size_t Ny, const std::string meshname)
{
    std::ofstream OUT(meshname.c_str());
    OUT << "ParametricMesh 1.0" << std::endl
        << Nx << "\t" << Ny << std::endl;
    for (size_t iy = 0; iy <= Ny; ++iy)
        for (size_t ix = 0; ix <= Nx; ++ix)
            OUT << Lx * ix / Nx << "\t" << Ly * iy / Ny << std::endl;
    OUT.close();
}

int main()
{
    Nextsim::ParametricMesh smesh(0); // 0 means no output

    for (int it = 0; it < 3; ++it) {
        create_rectanglemesh(ProblemConfig::Lx, ProblemConfig::Ly, ProblemConfig::Nx, ProblemConfig::Ny, "tmp.smesh");
        smesh.readmesh("tmp.smesh");

        Test<1> test1(smesh);
        test1.run();
        Test<3> test3(smesh);
        test3.run();
        Test<6> test6(smesh);
        test6.run();

        ProblemConfig::Nx *= 2;
        ProblemConfig::Ny *= 2;
    }
}
