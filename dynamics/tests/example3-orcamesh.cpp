/*!
 * @file eexample1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "DGTransport.hpp"
#include "dgLimiters.hpp"
#include "dgVisu.hpp"

#include "stopwatch.hpp"
#include "testtools.hpp"
#include <cassert>
#include <chrono>
#include <filesystem> // only for automatic creation of output directory
#include <iostream>
#include <vector>

bool WRITE_VTK = true; //!< set to true for vtk output
int WRITE_EVERY = 5;

double TOL = 1.e-10; //!< tolerance for checking test results

/*!
 *  Description of the test case
 *
 * Transport an initial smooth bump diagonally through the domain. 
 * The mesh is the NH-mesh projected to 2d (by setting z-coordinate to zero)
 * All calculations in Cartesian coordinates
 *
 */

namespace ProblemConfig {
  const double Lx = 0.5e6;
  const double Ly = 1.0e6;
  size_t Nx = -1;
  size_t Ny = -1;
  const double T  = 2.0 * M_PI;
  size_t NT = -1;
}


//! The initial solution
class SmoothBump : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
        double X = x / ProblemConfig::Lx;
        double Y = y / ProblemConfig::Ly;
        double r = (pow(0.5 * (X - 0.25), 2.0) + pow((Y - 0.25), 2.0)) / 0.02;
        if (r < 1)
            return exp(-1.0 / (1.0 - r));
        else
            return 0.0;
    }
};
//! F=0 for computing the error
class ZeroFunction : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
      return 0.0;
    }
};

// Velocity
class InitialVX : public Nextsim::Interpolations::Function { // (0.5,0.2) m/s

public:
    double operator()(double x, double y) const
    {
      return -(y-400000.0);
    }
};
class InitialVY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return (x-100000.0);
    }
};

//////////////////////////////////////////////////

template <int DG>
class Test {

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
      dt = ProblemConfig::T/ProblemConfig::NT;// 1000 sek.
      //! Init Vectors
      phi.resize_by_mesh(smesh);
    }

    double run()
    {

        //! Compose name of output directory and create it
        std::string resultsdir = "Example3_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);
        std::filesystem::create_directory(resultsdir);

	// init the test case, in particular resize vectors
        init();


        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, SmoothBump());

        // velocity field
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);

        if (WRITE_VTK) {
            Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", 0, phi, smesh);
        }
        std::cout << DG << "\t" << ProblemConfig::NT << "\t" << smesh.nx << "\t" << std::flush;

        //! time loop
        for (size_t iter = 1; iter <= ProblemConfig::NT; ++iter) {

            dgtransport.reinitnormalvelocity();
            dgtransport.step(dt, phi); // performs one time step with the 2nd or 3rd Order Heun scheme
            if (WRITE_VTK)
                if (iter % (ProblemConfig::NT / writestep) == 0)
		  {
                    Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", iter / (ProblemConfig::NT / writestep), phi, smesh);
		  }
        }
        // integral over the solution
        return Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, SmoothBump());
    }
};

template <int DG>
void run(double distort = 0.0)
{
    Nextsim::ParametricMesh smesh(Nextsim::CARTESIAN); // 0 means no output

#define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2 \
                                                         : -1))

    // Read the mesh
    smesh.readmesh("example3-25km_NH.smesh");
    ProblemConfig::Nx = smesh.nx;
    ProblemConfig::Ny = smesh.ny;
    ProblemConfig::NT = 500*(2*DG2DEG(DG)+1);
    // Transform mesh coordinates to Cartesian coordinates, then project to 2d
    double R = 6371000.0;  
    for (size_t i = 0; i<smesh.nnodes;++i)
      {
	const double x = R * cos(M_PI/180.0*smesh.vertices(i,1)) * cos(M_PI/180.0*smesh.vertices(i,0));
	const double y = R * cos(M_PI/180.0*smesh.vertices(i,1)) * sin(M_PI/180.0*smesh.vertices(i,0));
	smesh.vertices(i,0) =0.25*(-sin(M_PI/4.)*x+cos(M_PI/4.0)*y       );
	smesh.vertices(i,1) =0.25*( cos(M_PI/4.)*x+sin(M_PI/4.0)*y + 2.e6);
      }
    
    Test<DG> test(smesh);
    double integral = test.run();

    /*!
     * Exact values taken on 14.12.2022
     * We cannot expect convergence as the initial pattern
     * will cross some coastlines while rotating.
     */
    double exact[3]={4.9779942167207122e+08,1.3467162055206916e+06,1.4092155385515299e+06};

    std::cout << std::setprecision(4) << std::scientific;
    double error = exact[DG2DEG(DG)]-integral;
    double tol   = 1.e-8;
    bool passed = (fabs(error)/fabs(integral)<tol);
    std::cout << integral << "\t" << error << "\t"
	      << passed << std::endl;
}

int main()
{
  std::cout << "DG\tN\tNX\tIntegral\tError\t\tPassed (1)" << std::endl;
  run<1>();
  run<3>();
  run<6>();
  
    return 0;
}
