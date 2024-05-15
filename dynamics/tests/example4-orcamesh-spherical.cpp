/*!
 * @file eexample1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "Tools.hpp"
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

double TOL = 1.e-10; //!< tolerance for checking test results

/*!
 *  Description of the test case
 *
 * Comp. in spherical
 *
 */



namespace ProblemConfig {
  double T  = 2.0 * M_PI * Nextsim::EarthRadius;       // time horizon
  size_t NT = -1;         // no. time steps

  std::map<std::pair<size_t,size_t>, double> exact;
  
}


//! The initial solution
class SmoothBump : public Nextsim::Interpolations::Function {

public:
    double operator()(double lon, double lat) const
  {
    double rx = fabs(lon+2.0/3.0 * M_PI)/ (M_PI/6.0);
    double ry = fabs(lat) / (M_PI / 5.0);

    if (rx < 1 && ry < 1)
      return exp(-0.25 / (1.0 - rx*rx)) * exp(-0.25 / (1.0 - ry*ry))*exp(0.25*0.25);
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
      return cos(y); // 1 m/s at equator
    }
};
class InitialVY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return 0.0;
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
    Nextsim::DGTransport<DG> sphericaltransport;


    //! Velocity Field
    InitialVX VX;
    InitialVY VY;

    size_t writestep; //! write out n step in total (for debugging only)

public:
    Test(const Nextsim::ParametricMesh& mesh)
        : smesh(mesh)
        , sphericaltransport(smesh)
        , writestep(40)
    {
        //! Set time stepping scheme. 2nd order for dg0 and dg1, 3rd order dG2
        if (DG < 3)
            sphericaltransport.settimesteppingscheme("rk2");
        else
            sphericaltransport.settimesteppingscheme("rk3");
    }

    Test() { }

    void init()
    {
      dt = ProblemConfig::T / ProblemConfig::NT ;
      //! Init Vectors
      phi.resize_by_mesh(smesh);
    }

    double run()
    {

        //! Compose name of output directory and create it
        std::string resultsdir = "Example4_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);
        std::filesystem::create_directory(resultsdir);

	// init the test case, in particular resize vectors
        init();


        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, SmoothBump());
	
        // velocity field
        Nextsim::Interpolations::Function2DG(smesh, sphericaltransport.GetVx(), VX);
        Nextsim::Interpolations::Function2DG(smesh, sphericaltransport.GetVy(), VY);

	
        std::cout << DG << "\t" << ProblemConfig::NT << "\t" << smesh.nx << "\t" << std::flush;

	if (WRITE_VTK)
	  {
	    Nextsim::VTK::write_dg<1>(resultsdir + "/landmask", 0, Nextsim::Tools::Landmask(smesh), smesh);
	    Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", 0, phi, smesh);
	    Nextsim::VTK::write_dg<DG>(resultsdir + "/vx", 0, sphericaltransport.GetVx(), smesh);
	    Nextsim::VTK::write_dg<DG>(resultsdir + "/vy", 0, sphericaltransport.GetVy(), smesh);
	  }

	
        //! time loop
        for (size_t iter = 1; iter <= ProblemConfig::NT; ++iter) {

            sphericaltransport.reinitnormalvelocity();
            sphericaltransport.step(dt, phi); // performs one time step with the 2nd or 3rd Order Heun scheme
            if (WRITE_VTK)
                if (iter % (ProblemConfig::NT / writestep) == 0)
		  Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", iter / (ProblemConfig::NT / writestep), phi, smesh);
        }
        // integral over the solution
        return Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, SmoothBump());
    }
};


//! Creates a rectangular mesh of the whole earth. Ice for  poles * Ny < iy < (1-poles) * Ny
void create_rectanglemesh(const std::string meshname, size_t Nx, size_t Ny, double poles) 
{
    std::ofstream OUT(meshname.c_str());
    OUT << "ParametricMesh 2.0" << std::endl
        << Nx << "\t" << Ny << std::endl;
    for (size_t iy = 0; iy <= Ny; ++iy) // lat 
      for (size_t ix = 0; ix <= Nx; ++ix) // lon
	OUT << -M_PI+2.0 * M_PI*ix/Nx << "\t" << -M_PI/2.0+M_PI*iy/Ny << std::endl;

    OUT << "landmask " << Nx * Ny << std::endl; // no ice on poles :-)
    size_t y0 = poles*Ny;       // first element that is ice
    size_t y1 = (1.0-poles)*Ny; // last element that is ice
    for (size_t iy = 0; iy < Ny; ++iy) // lat 
      for (size_t ix = 0; ix < Nx; ++ix) // lon
	if ( (y0 <= iy) && (iy <= y1) )
	  OUT << 1 << std::endl;
	else
	  OUT << 0 << std::endl;

    // Dirichlet boundary along y0 (bottom) / y1 (top)
    OUT << "dirichlet " << 2*Nx << std::endl; // horizontal
    for (size_t i=0;i<Nx;++i)
      OUT << y0*Nx + i << "\t" << 0 << std::endl; // lower
    for (size_t i=0;i<Nx;++i)
      OUT << y1*Nx + i << "\t" << 2 << std::endl; // upper
    
    OUT << "periodic 1" << std::endl; // Periodic Y-term [1] left/right
    OUT << Ny << std::endl;
    for (size_t i=0;i<Ny;++i)
      OUT << (i+1)*Nx-1 << "\t" << i*Nx << "\t1" << std::endl;

    OUT.close();
}

template <int DG>
void run(size_t N)
{
  Nextsim::ParametricMesh smesh(Nextsim::SPHERICAL); // true = spherical 0 means no output

#define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2 \
                                                         : -1))

    // Read the mesh
  create_rectanglemesh("example4.smesh", 2*N, N, 0.2);
  smesh.readmesh("example4.smesh");
  // time step size
  ProblemConfig::NT = N*5*(2*DG2DEG(DG)+1);
  //  smesh.RotatePoleToGreenland();
    //    smesh.RotatePoleFromGreenland();

    Test<DG> test(smesh);
    double integral = test.run();

    std::cout << std::setprecision(4) << std::scientific;
    double tol   = 1.e-8;
    bool passed = (fabs(integral - ProblemConfig::exact[{DG,N}])/fabs(integral)<tol);
    std::cout << integral << "\t" << ProblemConfig::exact[{DG,N}] << "\t"
	      << passed << std::endl;
}

int main()
{
  // meshes WITHOUT rotation to Greenland. Values taken 15.01.2023 (radians & m/s) 
  ProblemConfig::exact[{1,32} ] = 1.2164748581905749e-01;
  ProblemConfig::exact[{3,32} ] = 5.6513958051583706e-03;
  ProblemConfig::exact[{6,32} ] = 9.5899046465648826e-04;
  ProblemConfig::exact[{1,64} ] = 8.1798148796320413e-02;
  ProblemConfig::exact[{3,64} ] = 1.7539002988624004e-03;
  ProblemConfig::exact[{6,64} ] = 1.6410486621918751e-04;
  ProblemConfig::exact[{1,128}] = 4.7832176725503828e-02;
  ProblemConfig::exact[{3,128}] = 4.4521911163719342e-04;
  ProblemConfig::exact[{6,128}] = 2.1342844522123369e-05;

  

  std::cout << "DG\tNT\tNX\tError\t\tReference\tPassed (1)" << std::endl;
  
  // for (size_t n : {32,64})
  for (size_t n : {32,64,128}) // Reference values for N=128 are given. But computations are lengthy.
    {
      run<1>(n);
      run<3>(n);
      run<6>(n);
    }
  
    return 0;
}
