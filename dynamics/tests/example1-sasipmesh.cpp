/*!
 * @file eexample1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "SphericalTransport.hpp"
#include "dgLimiters.hpp"
#include "dgVisu.hpp"
#include "Tools.hpp"


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
  size_t Nx = 24;
  const size_t Nx0 = 24;
  size_t Ny = 26;
  const size_t Ny0 = 26;
  size_t NT = 200;
  const size_t NT0 = 200;
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

    //! Compose name of output directory and create it
    std::string resultsdir = "Example1_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);
    std::filesystem::create_directory(resultsdir);

    // init the test case, in particular resize vectors
    init();

    // initial density
    Nextsim::Interpolations::Function2DG(smesh, phi, SmoothBump());

    double initialaverage = Nextsim::Tools::MeanValue(smesh,phi);

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
	  Nextsim::VTK::write_dg<DG>(resultsdir + "/dg", iter / (ProblemConfig::NT / writestep), phi, smesh);
    }

    initialaverage = (initialaverage - Nextsim::Tools::MeanValue(smesh,phi))/initialaverage; // compute error in mass
    std::cout << initialaverage << "\t";

    // integrate the error
    return sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, phi, SmoothBump())) / ProblemConfig::Lx;
  }
};

//////////////////////////////////////////////////

void create_rectanglemesh(const double Lx, const double Ly, const size_t Nx, const size_t Ny, const std::string meshname, double distort = 0.0)
{
  std::ofstream OUT(meshname.c_str());
  OUT << "ParametricMesh 2.0" << std::endl
      << Nx << "\t" << Ny << std::endl;
  for (size_t iy = 0; iy <= Ny; ++iy)
    for (size_t ix = 0; ix <= Nx; ++ix)
      OUT << Lx * ix / Nx + Lx * distort * sin(M_PI * ix / Nx * 3.0) * sin(M_PI * iy / Ny) << "\t"
	  << Ly * iy / Ny + Ly * distort * sin(M_PI * iy / Ny * 2.0) * sin(M_PI * ix / Nx * 2.0) << std::endl;

  // landmask
  OUT << "landmask 0" << std::endl;

  // dirichlet info
  OUT << "dirichlet " << 2*Nx + 2*Ny << std::endl;
  for (size_t i = 0; i < Nx; ++i)
    OUT << i << "\t" << 0 << std::endl; // lower

  for (size_t i = 0; i < Ny; ++i)
    OUT << i * Nx + Nx - 1 << "\t" << 1 << std::endl; // right
    
  for (size_t i = 0; i < Nx; ++i)
    OUT << Nx * (Ny - 1) + i << "\t" << 2 << std::endl; // upper
    
  for (size_t i = 0; i < Ny; ++i)
    OUT << i * Nx << "\t" << 3 << std::endl; // left

  OUT << "periodic 0" << std::endl;
  OUT.close();
}

template <int DG>
void run(double distort, const std::array<std::array<double, 4>, 3>& exact)
{
  Nextsim::ParametricMesh smesh; // 0 means no output

#define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2	\
				   : -1))

  for (int it = 0; it < 4; ++it) {
    ProblemConfig::Nx = ProblemConfig::Nx0 * (1 << it);
    ProblemConfig::Ny = ProblemConfig::Ny0 * (1 << it);
    ProblemConfig::NT = ProblemConfig::NT0 * (DG2DEG(DG) + 1) * (DG2DEG(DG) + 1) * (1 << it);

    create_rectanglemesh(ProblemConfig::Lx, ProblemConfig::Ly, ProblemConfig::Nx, ProblemConfig::Ny, "tmp1.smesh", distort);
    smesh.readmesh("tmp1.smesh");

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
  std::cout << std::setprecision(3) << std::scientific;
      
  std::array<std::array<double, 4>, 3> exact= // Exact values taken 14.12.2022 (with rk1/rk2/rk3)
    {
      std::array<double,4>({
	  4.8211028295934204e-02,
	  4.3147199226169552e-02,
	  3.5825194395857747e-02,
	  2.6956405683141109e-02}),
      std::array<double,4>({
	  1.3469076224946228e-02,
	  5.8684338697161317e-03,
	  2.2778998006559789e-03,
	  7.3625764726383837e-04}),
      std::array<double,4>({
	  4.0419818925870685e-03,
	  1.2394804214472853e-03,
	  2.8473714699424552e-04,
	  4.5484079480627463e-05})
    };

  std::cout << std::endl
	    << "DG\tNT\tNX\tmass loss\terror\t\texact\t\tpassed" << std::endl; 
  run<1>(0.0,exact);
  run<3>(0.0,exact);
  run<6>(0.0,exact);

  
  exact=
    {
      std::array<double,4>({
	  4.8651365310526981e-02,
	  4.3965656991140793e-02,
	  3.7061326881950421e-02,
	  2.8347869256577524e-02}),
      std::array<double,4>({
	  1.5618850741989440e-02,
	  6.7948174536363922e-03,
	  2.7387757081918577e-03,
	  9.2152265548481316e-04}),
      std::array<double,4>({
	  4.9194576782567610e-03,
	  1.5855619683891836e-03,
	  3.8872699850545491e-04,
	  6.7179605760583785e-05})
    };
  
  std::cout << std::endl << "Distorted mesh" << std::endl;
  std::cout << "DG\tNT\tNX\tmass loss\terror\t\texact\t\tpassed" << std::endl;
  run<1>(0.05,exact);
  run<3>(0.05,exact);
  run<6>(0.05,exact);


  return 0;
}
