/*!
 * @file test1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "DGTransport.hpp"
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
      Nextsim::VTK::write_dg<DG>(resultsdir + "/vx", 0, dgtransport.GetVx(), smesh);
      Nextsim::VTK::write_dg<DG>(resultsdir + "/vy", 0, dgtransport.GetVy(), smesh);
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
  Nextsim::ParametricMesh smesh(Nextsim::CARTESIAN); // 0 means no output

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
  std::cout << std::setprecision(4) << std::scientific;

  std::array<std::array<double, 4>, 3> exact= // Exact values taken 14.12.2022 (with rk1/rk2/rk3)
    {
      std::array<double,4>({
	  4.8314783396276907e-02,
	  4.3121845060990968e-02,
	  3.5805596896499432e-02,
	  2.6943655410656617e-02}),
      std::array<double,4>({
	  1.3427300428974142e-02,
	  5.8574917975990539e-03,
	  2.2756828198281097e-03,
	  7.3564123952743888e-04}),
      std::array<double,4>({
	  4.0274885631835103e-03,
	  1.2400192553115746e-03,
	  2.8460172952094940e-04,
	  4.5459584746490405e-05})
    };

  std::cout << std::endl
	    << "DG\tNT\tNX\tmass loss\terror\t\texact\t\tpassed" << std::endl;
  run<1>(0.0,exact);
  run<3>(0.0,exact);
  run<6>(0.0,exact);

  exact=
    {
      std::array<double,4>({
	  4.8471012332113990e-02,
	  4.3923292993805568e-02,
	  3.7043103484653367e-02,
	  2.8337995655005246e-02}),/*!
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
    #include "Tools.hpp"


    #include "stopwatch.hpp"
    #include "testtools.hpp"
    #include <cassert>
    #include <chrono>
    #include <filesystem> // only for automatic creation of output directory
    #include <iostream>
    #include <vector>


    bool WRITE_VTK = false; //!< set to true for vtk output
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

        double initialaverage = Nextsim::Tools::MeanValue(smesh,phi);

        // velocity field
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
        Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);

        // Write out initial solution
        std::string resultsdir = "Example1_" + std::to_string(DG) + "_" + std::to_string(smesh.nx);
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

    void create_rectanglemesh(Nextsim::ParametricMesh& smesh, const double Lx, const double Ly, const size_t Nx, const size_t Ny, double distort = 0.0)
    {
      smesh.reset(); // reset mesh and clear all variables

      // set number of elements and nodes
      smesh.nx = Nx;
      smesh.ny = Ny;
      smesh.nnodes = (Nx+1)*(Ny+1);
      smesh.nelements = Nx * Ny;

      // set vertices
      smesh.vertices.resize(smesh.nnodes,2);
      size_t ii = 0;
      for (size_t iy = 0; iy <= Ny; ++iy)
        for (size_t ix = 0; ix <= Nx; ++ix, ++ii)
          {
    	smesh.vertices(ii,0) =  Lx * ix / Nx + Lx * distort * sin(M_PI * ix / Nx * 3.0) * sin(M_PI * iy / Ny);
    	smesh.vertices(ii,1) =  Ly * iy / Ny + Ly * distort * sin(M_PI * iy / Ny * 2.0) * sin(M_PI * ix / Nx * 2.0);
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

    template <int DG>
    void run(double distort, const std::array<std::array<double, 4>, 3>& exact)
    {
      Nextsim::ParametricMesh smesh(Nextsim::CARTESIAN); // 0 means no output

    #define DG2DEG(DG) (DG == 1 ? 0 : (DG == 3 ? 1 : DG == 6 ? 2	\
    				   : -1))

      for (int it = 0; it < 2; ++it) {

        // set number of elements in x- and y-direction
        ProblemConfig::Nx = ProblemConfig::Nx0 * (1 << it);
        ProblemConfig::Ny = ProblemConfig::Ny0 * (1 << it);
        // set number of time-steps such that the cfl condition is fulfilled
        ProblemConfig::NT = ProblemConfig::NT0 * (DG2DEG(DG) + 1) * (DG2DEG(DG) + 1) * (1 << it);

        // create the mesh
        create_rectanglemesh(smesh, ProblemConfig::Lx, ProblemConfig::Ly, ProblemConfig::Nx, ProblemConfig::Ny, distort);

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
      std::cout << std::setprecision(4) << std::scientific;

      std::array<std::array<double, 4>, 3> exact= // Exact values taken 14.12.2022 (with rk1/rk2/rk3)
        {
          std::array<double,4>({
    	  4.8314783396276907e-02,
    	  4.3121845060990968e-02,
    	  3.5805596896499432e-02,
    	  2.6943655410656617e-02}),
          std::array<double,4>({
    	  1.3427300428974142e-02,
    	  5.8574917975990539e-03,
    	  2.2756828198281097e-03,
    	  7.3564123952743888e-04}),
          std::array<double,4>({
    	  4.0274885631835103e-03,
    	  1.2400192553115746e-03,
    	  2.8460172952094940e-04,
    	  4.5459584746490405e-05})
        };

      std::cout << std::endl
    	    << "DG\tNT\tNX\tmass loss\terror\t\texact\t\tpassed" << std::endl;
      run<1>(0.0,exact);
      run<3>(0.0,exact);
      run<6>(0.0,exact);

      exact=
        {
          std::array<double,4>({
    	  4.8471012332113990e-02,
    	  4.3923292993805568e-02,
    	  3.7043103484653367e-02,
    	  2.8337995655005246e-02}),
          std::array<double,4>({
    	  1.5639053380182857e-02,
    	  6.7927021150567882e-03,
    	  2.7369920558767478e-03,
    	  9.2109978015705776e-04}),
          std::array<double,4>({
    	  4.9207622570231089e-03,
    	  1.5831521586435636e-03,
    	  3.8869573083901034e-04,
    	  6.7162994609144087e-05})
        };


      std::cout << std::endl << "Distorted mesh" << std::endl;
      std::cout << "DG\tNT\tNX\tmass loss\terror\t\texact\t\tpassed" << std::endl;
      run<1>(0.05,exact);
      run<3>(0.05,exact);
      run<6>(0.05,exact);


      return 0;
    }

      std::array<double,4>({
	  1.5639053380182857e-02,
	  6.7927021150567882e-03,
	  2.7369920558767478e-03,
	  9.2109978015705776e-04}),
      std::array<double,4>({
	  4.9207622570231089e-03,
	  1.5831521586435636e-03,
	  3.8869573083901034e-04,
	  6.7162994609144087e-05})
    };


  std::cout << std::endl << "Distorted mesh" << std::endl;
  std::cout << "DG\tNT\tNX\tmass loss\terror\t\texact\t\tpassed" << std::endl;
  run<1>(0.05,exact);
  run<3>(0.05,exact);
  run<6>(0.05,exact);


  return 0;
}
