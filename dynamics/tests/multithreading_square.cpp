#include <cassert>
#include <iostream>
#include <chrono>

#include "dginitial.hpp"
#include "dglimiters.hpp"
#include "mesh.hpp"
#include "timemesh.hpp"
#include "dgvisu.hpp"
#include "dgtransport.hpp"
#include "stopwatch.hpp"


namespace Nextsim
{
  extern Timer GlobalTimer;
}


bool WRITE_VTK = false;


//! Initially a smooth bump centered at (0.4,0.4)
//! This will be transported in a circle with (-y, x) for one complete revolution
class InitialPhi : virtual public Nextsim::InitialBase {
public:
  double smooth(double x) const // smooth transition from 0 to 1 on [0,1]
  {
    if (x<=0)
      return 0;
    if (x>=1.0)
      return 0;
    
    if (x<0.5)
      return 0.5*exp(-1./x)/exp(-2.0);
    else
      return 1.-0.5*exp(-1./(1.-x))/exp(-2.0);
  }
  

  double operator()(double x, double y) const
  {
    double r = sqrt(pow(x-0.4,2.0)+pow(y-0.4,2.0));
    if (r<0.1)
      return 1.0;
    if (r<0.3)
      return 1.0-smooth(5.0*(r-0.1));
    return 0.0;
  }
};


// Velocity
class InitialVX : virtual public Nextsim::InitialBase {
public:
  virtual double operator()(double x, double y) const { return y - 0.5; }
};
class InitialVY : virtual public Nextsim::InitialBase {
public:
  virtual double operator()(double x, double y) const { return 0.5 - x; }
};

//////////////////////////////////////////////////


template<int DGdegree>
class Test
{
  //! Mesh size (given as parameter to constructor)
  size_t N;
  
  //! Meshes
  Nextsim::Mesh mesh;
  Nextsim::TimeMesh timemesh;
  
  //! Velocity vectors and density
  Nextsim::CellVector<DGdegree> vx, vy, phi, finalphi;

  //! Transport main class
  Nextsim::DGTransport<DGdegree> dgtransport;

  //! Velocity Field
  InitialVX VX;
  InitialVY VY;

 public:

  Test(size_t n) : N(n), dgtransport(vx,vy)
  {
    dgtransport.settimesteppingscheme("rk3");
  }

  Test() 
  {
    std::cerr << "call Test(N). N is number of mesh elements per row" << std::endl;
  }

  void init()
  {
    //! Init Mesh
    mesh.BasicInit(N, N, 1.0 / N, 1.0 / N);

    //! Init Time Mesh
    double cfl = 0.1;
    double k   = cfl * std::min(mesh.hx, mesh.hy) / 1.0; // max-velocity is 1
    double tmax = 2.0*M_PI;

    int NT = (static_cast<int>((tmax / k + 1) /100 + 1) * 100); // No time steps dividable by 100
    timemesh.BasicInit(tmax, NT);
    
    //! Init Transport Scheme
    dgtransport.setmesh(mesh);
    dgtransport.settimemesh(timemesh);

    //! Init Vectors
    vx.resize_by_mesh(mesh);
    vy.resize_by_mesh(mesh);
    phi.resize_by_mesh(mesh);
  }
  
  void run()
  {
    Nextsim::GlobalTimer.reset();
    Nextsim::GlobalTimer.start("--> run");


    Nextsim::GlobalTimer.start("-- --> initial");
    // initial density
    Nextsim::L2ProjectInitial(mesh, phi, InitialPhi());

    if (WRITE_VTK)
      Nextsim::VTK::write_dg<DGdegree>("Results/dg",0,phi,mesh);

    //! Save initial solution for error control
    finalphi = phi;

    // set velocity vector. constant velocity field. No initialization required
    Nextsim::L2ProjectInitial(mesh, vx, VX);
    Nextsim::L2ProjectInitial(mesh, vy, VY);
    dgtransport.reinitvelocity();

    Nextsim::GlobalTimer.stop("-- --> initial");
    // time loop
    for (size_t iter = 1; iter <= timemesh.N; ++iter)
      {
	dgtransport.step(phi); // performs one time step with the 2nd Order Heun scheme
	if (WRITE_VTK)
	  if (iter % (timemesh.N/10)==0)
	    Nextsim::VTK::write_dg<DGdegree>("Results/dg",iter/(timemesh.N/10),phi,mesh);
      }

    Nextsim::GlobalTimer.stop("--> run");

    Nextsim::GlobalTimer.print();
  }

  bool check() const
  {
    std::cerr << "k " << 1./N << " dG " << DGdegree << "\t";

    //! Check that mass is ok.
    double mass      = phi.mass(mesh);
    std::cerr << std::setprecision(16)
	      << mass << "\t";

    Nextsim::CellVector<DGdegree> errorphi = phi; errorphi += -finalphi;
    if (WRITE_VTK)
      Nextsim::VTK::write_dg<DGdegree>("Results/error",N,errorphi,mesh);
    
    std::cerr << errorphi.norm() * sqrt(mesh.hx*mesh.hy) << std::endl;
    
    return true;
  }
  

};


//////////////////////////////////////////////////

int main()
{
  
  int ITER = 3;
  int N = 20;
  for (int n=0;n<ITER;++n)
    {
      Test<0> test0(N);
      test0.init();
      test0.run();
      if (!test0.check())
  	std::cerr << "TEST FAILED!" << std::endl;
      N*=2;
    }
  std::cerr << std::endl;
  
  N = 20;
  for (int n=0;n<ITER;++n)
    {
      Test<1> test1(N);
      test1.init();
      test1.run();
      if (!test1.check())
  	std::cerr << "TEST FAILED!" << std::endl;
      N*=2;
    }
  std::cerr << std::endl;
  
  N = 20;
  for (int n=0;n<ITER;++n)
    {
      Test<2> test2(N);
      test2.init();
      test2.run();
      if (!test2.check())
	std::cerr << "TEST FAILED!" << std::endl;
      N*=2;
    }
}
