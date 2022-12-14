/*!
 * @file eexample1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "SphericalTransport.hpp"
#include "dgLimit.hpp"
#include "dgVisu.hpp"
#include "Tools.hpp"

#include "stopwatch.hpp"
#include "testtools.hpp"
#include <cassert>
#include <chrono>
#include <filesystem> // only for automatic creation of output directory
#include <iostream>
#include <vector>

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


//! exact error values. [dg][it]. right order of convergence (1/2), (2) and (3) need at least 4 meshes
double exact_values[6][4] = {
    { 0.0482094, 0.0431458, 0.0358241, 0.0269558 },
    { 0.0134817, 0.00587354, 0.00227998, 0.000737005 },
    { 0.00410855, 0.00127682, 0.000304336, 5.55225e-05 },
    { 0.0482094, 0.0431458, 0.0358241, 0.0269558 },
    { 0.0134817, 0.00587354, 0.00227998, 0.000737005 },
    { 0.00410855, 0.00127682, 0.000304336, 5.55225e-05 }
};

bool WRITE_VTK = false; //!< set to true for vtk output

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
	  return exp(1.0) * exp(-1.0/(1.0-r));

	r = (pow((x) / 50000, 2.0) + pow((y-175000.0) / 50000.0, 2.0));
	if (r < 1.0)
	  {
	    if (y>175000)
	      return 1.0;
	    if (fabs(x)>15000.0)
	      return 1.0;
	    return 0.0;
	  }

	r = sqrt(pow((x - 175000.0) / 50000, 2.0) + pow((y) / 50000.0, 2.0));
	if (r<1.0)
	  return 1.0-r;

	r = sqrt(pow((x) / 50000, 2.0) + pow((y+175000.) / 50000.0, 2.0));
	if (r<1.0)
	  return 1.0;


	
	
	else return 0.0;
    }
};

// Velocity
class InitialVX : public Nextsim::Interpolations::Function { // (0.5,0.2) m/s

public:
    double operator()(double x, double y) const
    {
        return y * 2.0 * M_PI / ProblemConfig::R1;
    }
};
class InitialVY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
        return - x * 2.0 * M_PI / ProblemConfig::R1;
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
    Nextsim::SphericalTransport<DG> dgtransport;

    //! Velocity Field
    InitialVX VX;
    InitialVY VY;

    size_t writestep; //! write out n step in total (for debugging only)

public:
    Test(const Nextsim::ParametricMesh& mesh)
        : smesh(mesh)
        , dgtransport(smesh, Nextsim::CARTESIAN)
        , writestep(1)
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
        //! Compose name of output directory and create it
        std::string resultsdir = "Example2_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);
        std::filesystem::create_directory(resultsdir);

        // init the test case, in particular resize vectors
        init();

        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, Packman());
	Nextsim::LimitMax(phi,1.0);
	Nextsim::LimitMin(phi,0.0);

	double initialmass = Nextsim::Tools::MeanValue(smesh,phi);
        // velocity field
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);

	
        if (WRITE_VTK) {
            Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", 0, phi, smesh);
        }
        std::cout << DG << "\t" << ProblemConfig::NT << "\t" << smesh.nx << "\t" << std::flush;


        // time loop
        for (size_t iter = 1; iter <= ProblemConfig::NT; ++iter) {

            dgtransport.reinitnormalvelocity();
            dgtransport.step(dt, phi); // performs one time step with the 2nd or 3rd Order Heun scheme
	    Nextsim::LimitMax(phi,1.0);
	    Nextsim::LimitMin(phi,0.0);

            if (WRITE_VTK)
	      if (iter % (ProblemConfig::NT / writestep) == 0)
                    Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", iter / (ProblemConfig::NT / writestep), phi, smesh);
        }
	initialmass = (initialmass - Nextsim::Tools::MeanValue(smesh,phi)) / initialmass;
	
	std::cout << initialmass << "\t";
	// compute error
	return sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, Packman())) / ProblemConfig::R0;
    }
};

//////////////////////////////////////////////////

void create_rectanglemesh(const std::string meshname)
{
    std::ofstream OUT(meshname.c_str());
    OUT << "ParametricMesh 2.0" << std::endl
        << ProblemConfig::Nx << "\t" << ProblemConfig::Ny << std::endl;
    for (size_t iy = 0; iy <= ProblemConfig::Ny; ++iy)
        for (size_t ix = 0; ix <= ProblemConfig::Nx; ++ix)
	  {
	    double r = ProblemConfig::R0 + (ProblemConfig::R1 - ProblemConfig::R0) * iy / ProblemConfig::Ny;
	    double p = -2.0 * M_PI * ix/ProblemConfig::Nx;
	    OUT << r*cos(p) << "\t" << r*sin(p) << std::endl;
	  }

    OUT << "landmask 0" << std::endl; // all ice
    
    OUT << "dirichlet " << 2*ProblemConfig::Nx << std::endl; // horizontal
    for (size_t i=0;i<ProblemConfig::Nx;++i)
      OUT << i << "\t" << 0 << std::endl; // lower
    for (size_t i=0;i<ProblemConfig::Nx;++i)
      OUT << ProblemConfig::Nx*(ProblemConfig::Ny-1)+i << "\t" << 2 << std::endl; // upper
    
    OUT << "periodic 1" << std::endl; // Periodic Y-term [1] left/right
    OUT << ProblemConfig::Ny << std::endl;
    for (size_t i=0;i<ProblemConfig::Ny;++i)
      OUT << (i+1)*ProblemConfig::Nx-1 << "\t" << i*ProblemConfig::Nx << "\t1" << std::endl;

    OUT.close();
}

template <int DG>
void run(const std::array<std::array<double, 4>, 3>& exact)
{
    Nextsim::ParametricMesh smesh; 

#define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2 \
                                                         : -1))

    for (int it = 0;it<4;++it)
      {
	
	ProblemConfig::Nx = ProblemConfig::Nx0 * (1 << it);
	ProblemConfig::Ny = ProblemConfig::Ny0 * (1 << it);
	ProblemConfig::NT = ProblemConfig::NT0 * (DG2DEG(DG) + 1) * (DG2DEG(DG) + 1) * (1 << it);
	
	create_rectanglemesh("tmp2.smesh");
	smesh.readmesh("tmp2.smesh");
	
	
	Test<DG> test(smesh);
	double error = test.run();
	std::cout << error << "\t" << exact[DG2DEG(DG)][it];
	if (fabs(error - exact[DG2DEG(DG)][it]) / exact[DG2DEG(DG)][it] < 1.e-4)
	  std::cout << "\ttest passed" << std::endl;
	else
	  std::cout << "\ttest failed" << std::endl;

	 
      }
}

int main()
{
  std::cout << "DG\tNT\tNX\tmass\t\terror\t\texact\t\tpassed" << std::endl;
  std::cout << std::setprecision(4) << std::scientific;

  std::array<std::array<double, 4>, 3> exact= // Exact values taken 22.11.2022
    {
      std::array<double,4>({
	  1.1694670693410663e+00,
	  1.1285759937350424e+00,
	  1.0659725870386016e+00,
	  9.5109023920919167e-01
	}),
      std::array<double,4>({
	  1.0566915186500225e+00,
	  7.6902733263023504e-01,
	  5.1638380542567819e-01,
	  3.6015749087091753e-01
	}),
      std::array<double,4>({
	  6.5986370457896470e-01,
	  4.1595435088567217e-01,
	  3.0395405738986181e-01,
	  2.2232746207904747e-01
	})
    };
  
  run<1>(exact);
  run<3>(exact);
  run<6>(exact);
  return 0;
}
