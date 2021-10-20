#include <cassert>
#include <iostream>
#include <iomanip>

#include "dginitial.hpp"
#include "mesh.hpp"
#include "timemesh.hpp"
#include "dgvisu.hpp"
#include "dgtransport.hpp"

bool WRITE_VTK = false;

/*!
 * This test case tests the boundary  handling of the DG transport scheme
 * An initial density is first transported to the upper right corner, 
 * then back to the lower left and finally back to the origin.
 * All boundaries are involved and we check if the final solution has the correct 
 * behavior
 *
 * Domain is [0,2] x [0,1] with 2N * N elements
 */

class InitialVX : virtual public Nextsim::InitialBase {
  double time;
  
public:

  void settime(double t)
  {
    time = t;
  }
  
  virtual double operator()(double x, double y) const
  {
    if ((time<0.4) || (time>1.2)) return 2.;
    return -2.;
  }
};
class InitialVY : virtual public Nextsim::InitialBase {
  double time;
public:
  void settime(double t)
  {
    time = t;
  }
  
  virtual double operator()(double x, double y) const
  {
    if ((time<0.4) || (time>1.2)) return 1.;
    return -1.;
  }
};
class InitialPhi : virtual public Nextsim::InitialBase {
public:
    virtual double operator()(double x, double y) const
  {
    return exp(-50.0 * pow(x-1.0,2.0) - 50.0 * pow(y-0.5,2.0));
  }
};

template<int DGdegree>
class Test
{
  //! Meshes
  Nextsim::Mesh mesh;
  Nextsim::TimeMesh timemesh;
  
  //! Velocity vectors and density
  Nextsim::CellVector<DGdegree> vx, vy, phi;

  //! Transport main class
  Nextsim::DGTransport<DGdegree> dgtransport;

  //! Velocity Field
  InitialVX VX;
  InitialVY VY;

 public:

  Test() : dgtransport(vx,vy)
  {
    dgtransport.settimesteppingscheme("rk2");
  }

  //! Returns the reference values for N=50 obtained on October 16, 2021
  double reference() const
  {
    assert(mesh.nx == 100);
    assert(mesh.ny == 50);
    
    if (DGdegree == 0)
      return 0.0209529644049705;
    else if (DGdegree == 1)
      return 0.04068581235936104;
    else if (DGdegree ==2)
      return 0.04075915956821777;
    abort();
    
  }
  
  void init()
  {
    //! Init Mesh
    size_t N = 50;
    mesh.BasicInit(2 * N, N, 1.0 / N, 1.0 / N);

    //! Init Time Mesh
    double cfl = 0.1; // 0.1 is the value used to get the reference values above
    
    double k   = cfl * std::min(mesh.hx, mesh.hy) / 2.0; // max-velocity is 1
    double tmax = 1.6;

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
    // initial density
    Nextsim::L2ProjectInitial(mesh, phi, InitialPhi());

    // time loop
    for (size_t iter = 1; iter <= timemesh.N; ++iter)
      {
	// set velocity vector
	VX.settime(iter * timemesh.k);
	VY.settime(iter * timemesh.k);
	Nextsim::L2ProjectInitial(mesh, vx, VX);
	Nextsim::L2ProjectInitial(mesh, vy, VY);

	dgtransport.reinitvelocity();// sets the current velocity and averages it to edges
	
	dgtransport.step(phi); // performs one time step with the 2nd Order Heun scheme
	if (WRITE_VTK)
	  if (iter % (timemesh.N/10)==0)
	    Nextsim::VTK::write_dg<DGdegree>("Results/dg",iter/(timemesh.N/10),phi,mesh);

      }
  }

  bool check() const
  {
    // integral over the [0.8,1.2] x [0.4,0.6]
    double exactmass = 0.4094292816e-1;
    double refmass   = reference();
    double mass      = phi.mass(mesh);
    double masserror = fabs(exactmass-mass);


    
    std::cerr << "Mass [Exact / Reference / Numerical / Error]\t"
	      << std::setprecision(8)
	      << exactmass << "\t"
      	      << refmass << "\t"
	      << mass << "\t"
	      << masserror << "\t"
	      << std::endl;
    return (fabs(mass - refmass) < 1.e-8);
  }
  

};


int main()
{
  Test<0> test0;
  test0.init();
  test0.run();
  if (!test0.check())
    std::cerr << "TEST FAILED!" << std::endl;
  
  Test<1> test1;
  test1.init();
  test1.run();
  if (!test1.check())
    std::cerr << "TEST FAILED!" << std::endl;


  Test<2> test2;
  test2.init();
  test2.run();
  if (!test2.check())
    std::cerr << "TEST FAILED!" << std::endl;
  

}
