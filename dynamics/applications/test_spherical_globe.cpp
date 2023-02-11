/*!
 * @file benchmark_cartesian_globe.cpp
 * @date 24 July 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

/*!
 *
 * Ice dynamics test case the globe in spherical coordinates
 * of the globe. 

 * Ice on [0,2 Pi] * [- 3/8 Pi, 3/8 Pi]
 * Periodic at x=0 and x = 2pi, Dirichlet at y = +/- 10000
 * 
 * Wind is perturbation of (10,0)
 * Ocean (1,0) * (1.0 - (x/ (3/8 Pi)^2) and zero on land
 *
 * Full ice cover, A=1 and H=0.5 
 */

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "DGTransport.hpp"

#include "Tools.hpp"
#include "cgParametricMomentum.hpp"
#include "cgVector.hpp"
#include "dgInitial.hpp"
#include "dgLimit.hpp"
#include "dgVisu.hpp"
#include "mevp.hpp"
#include "stopwatch.hpp"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>
#include "VectorManipulations.hpp"
#include <map>

bool WRITE_VTK = true;

/*!
 * Exact values for the test case taken on Feb 11, 2023
 * See below in main for reference values
 */

std::map<std::array<size_t,2>, std::array<double,5> > test1;
std::map<std::array<size_t,2>, std::array<double,2> > test2;

std::string check(double a, double b)
{
  if (fabs(a-b)<1.e-12)
    return " ok ";
  else
    return " EE ";
}

namespace Nextsim {
extern Timer GlobalTimer;
}

inline constexpr double SQR(double x) { return (x * x); }

namespace ReferenceScale {
  constexpr double R1 = Nextsim::EarthRadius * 0.8; // For z=R1, Dirichlet zero
}




// The exact solution
class VX : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return SQR(Nextsim::EarthRadius) * sin(x)*cos(y)*cos(y) * (-cos(y)*cos(y)+0.36);
    }
};
class VY : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return SQR(Nextsim::EarthRadius) * sin(x)*cos(y)*cos(y)*cos(y)*cos(x) * (cos(y)*cos(y)-0.36);
    }
};
// The exact strain rate tensor
class exE11 : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return Nextsim::EarthRadius * cos(x)*cos(y) *
	( -pow(cos(y),3.0)*sin(x)*sin(y)+0.36*cos(y)*sin(x)*sin(y)-cos(y)*cos(y)+0.36);
    }
};
class exE12 : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return Nextsim::EarthRadius * cos(y) *
	((cos(x)*cos(x)-0.5)*pow(cos(y),3.0)+1.5*cos(y)*cos(y)*sin(x)*sin(y)+(0.18-0.36*cos(x)*cos(x))*cos(y)-0.18*sin(y)*sin(x));
    }
};
class exE22 : public Nextsim::Interpolations::Function {
public:
    double operator()(double x, double y) const
    {
      return Nextsim::EarthRadius * cos(y)*cos(y)*sin(x)*sin(y)*cos(x) *
	(1.08-5.0*cos(y)*cos(y));
    }
};

// The laplace
class LX : public Nextsim::Interpolations::Function {
public:
  double operator()(double X, double Y) const
  {
    const double t1 = cos(Y);
    const double t2 = t1 * t1;
    const double t3 = t2 * t2;
    const double t4 = sin(X);
    const double t7 = cos(X);
    const double t8 = t7 * t7;
    const double t11 = sin(Y);
    const double t24 = 0.9e1 * t4 * t3 + t2 * t1 * t11 * (-0.400e3 * t8 + 0.200e3) / 0.50e2 - 0.361e3 / 0.50e2 * t4 * t2 + t1 * t11 * (0.108e3 * t8 - 0.54e2) / 0.50e2 + 0.9e1 / 0.50e2 * t4;
    return t24;
  }
};
class LY : public Nextsim::Interpolations::Function {
public:
  double operator()(double X, double Y) const
  {
    const double t1 = cos(X);
    const double t2 = cos(Y);
    const double t3 = t2 * t2;
    const double t4 = t3 * t3;
    const double t6 = sin(X);
    const double t12 = sin(Y);
    const double t20 = (0.1298e4 * t6 * t3 * t2 - 0.1450e4 * t6 * t4 * t2 + 0.25e2 * t12 * t3 - 0.108e3 * t6 * t2 + 0.9e1 * t12) * t1 / 0.50e2;
    return t20;
  }
};
//////////////////////////////////////////////////


//! Creates a rectangular mesh of (nearly) whole earth. Ice for  poles * Ny < iy < (1-poles) * Ny
void create_mesh(Nextsim::ParametricMesh& smesh, size_t Nx, size_t Ny) 
{
  smesh.statuslog = -1;
  smesh.CoordinateSystem = Nextsim::SPHERICAL;
  
  smesh.nx = Nx;
  smesh.ny = Ny;
  smesh.nelements = Nx*Ny;
  smesh.nnodes    = (Nx+1)*(Ny+1);
  smesh.vertices.resize(smesh.nnodes,2);
  
  // z coordinate between -0.8*R and 0.8*R
  const double thetamax = asin(0.8);
  
  for (size_t iy = 0; iy <= Ny; ++iy) // lat 
    for (size_t ix = 0; ix <= Nx; ++ix) // lon
      {
	smesh.vertices(iy*(Nx+1)+ix,0) = -M_PI+2.0 * M_PI*ix/Nx;
	smesh.vertices(iy*(Nx+1)+ix,1) = -thetamax+2.0*thetamax*iy/Ny;
      }

  // ice everywhere
  smesh.landmask.resize(Nx*Ny);
  for (size_t i=0;i<Nx*Ny;++i)
    smesh.landmask[i]=1;

  // dirichlet boundary
  for (auto &it :  smesh.dirichlet)
    it.clear();
  for (size_t i=0;i<Nx;++i)
    {
      smesh.dirichlet[0].push_back(i);
      smesh.dirichlet[2].push_back((Ny-1)*Nx + i);
    }

  // periodic boundary
  smesh.periodic.clear();
  smesh.periodic.resize(1); // 1 segments
  for (size_t i=0;i<Ny;++i)
    smesh.periodic[0].push_back(std::array<size_t,4> ({1, (i+1)*Nx-1, i*Nx, i*(Nx+1)}));
}

template <int CG, int DGadvection>
void run_benchmark_strain(const size_t N)
{
    //! Define the spatial mesh
  Nextsim::ParametricMesh smesh(Nextsim::SPHERICAL);
  create_mesh(smesh, 2*N, N);

    //! Compose name of output directory and create it
    std::string resultsdir = "TestGlobe_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(N);
    std::filesystem::create_directory(resultsdir);
    
    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh, Nextsim::SPHERICAL);


    std::cout << CG << "\t" << N << "\t";

    ////////////////////////////////////////////////// Set Velocity
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetVx(), VX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetVy(), VY());
    momentum.DirichletZero(momentum.GetVx());
    momentum.DirichletZero(momentum.GetVy());
    Nextsim::VectorManipulations::CGAveragePeriodic(smesh, momentum.GetVx());
    Nextsim::VectorManipulations::CGAveragePeriodic(smesh, momentum.GetVy());

    // (1) Write out initial velocity
    Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", 0, momentum.GetVx(), momentum.GetVy(), smesh);

    Nextsim::DGVector<CG2DGSTRESS(CG)> E11(smesh),E12(smesh),E22(smesh);
    Nextsim::Interpolations::Function2DG(smesh, E11, exE11());
    Nextsim::Interpolations::Function2DG(smesh, E12, exE12());
    Nextsim::Interpolations::Function2DG(smesh, E22, exE22());
 
    
    // (2) Project Velocity to DG strain rate tensor
    momentum.ProjectCGVelocityToDGStrain();
    // compute error w.r.t. exact strain rate tensor given analytically (scale by 1/R)
    double e11 = sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, momentum.GetE11(), exE11()))/Nextsim::EarthRadius; 
    double e12 = sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, momentum.GetE12(), exE12()))/Nextsim::EarthRadius;
    double e22 = sqrt(Nextsim::Interpolations::L2ErrorFunctionDG(smesh, momentum.GetE22(), exE22()))/Nextsim::EarthRadius;

    
    

    // (3) Compute the divergence of the strain rate tensor (sym vector laplace beltrami)
    momentum.GetVx().zero();
    momentum.GetVy().zero();
    momentum.GetS11() = momentum.GetE11();
    momentum.GetS12() = momentum.GetE12();
    momentum.GetS22() = momentum.GetE22();
    Nextsim::CGVector<CG> tmpx(smesh), tmpy(smesh);
    momentum.DivergenceOfStress(1.0, tmpx, tmpy);
#pragma omp parallel for
    for (size_t i=0;i<momentum.GetVx().rows();++i)
      {
	momentum.GetVx()(i) = tmpx(i) / momentum.pmap.lumpedcgmass(i);
	momentum.GetVy()(i) = tmpy(i) / momentum.pmap.lumpedcgmass(i);
      }

    // For comparison get the exact laplacian
    Nextsim::CGVector<CG> laplaceX(smesh),laplaceY(smesh);
    Nextsim::Interpolations::Function2CG(smesh, laplaceX, LX());
    Nextsim::Interpolations::Function2CG(smesh, laplaceY, LY());
    Nextsim::VectorManipulations::CGAveragePeriodic(smesh,laplaceX);
    Nextsim::VectorManipulations::CGAveragePeriodic(smesh,laplaceY);
    // to compute error exclude boundary layer by adding one landmask-layer
    for (size_t ix=0;ix<smesh.nx;++ix)
      {
	smesh.landmask[                        ix]=0;
	smesh.landmask[smesh.nx*(smesh.ny-1) + ix]=0;
      }
    double l1 = sqrt(Nextsim::Interpolations::L2ErrorFunctionCG(smesh, momentum.GetVx(), LX()));
    double l2 = sqrt(Nextsim::Interpolations::L2ErrorFunctionCG(smesh, momentum.GetVy(), LY()));

    // write out discrete laplace and difference to exact one
    // first, set to zero on new landmasks
    size_t lastnode = (CG*smesh.nx+1)*(CG*smesh.ny+1)-1;
    for (size_t ix=0;ix<CG*smesh.nx+1;++ix)
      for (size_t iy=0;iy<CG+1;++iy)
	{
	  momentum.GetVx()(ix + (CG*smesh.nx+1) * iy)=0.0;
	  momentum.GetVy()(ix + (CG*smesh.nx+1) * iy)=0.0;
	  laplaceX        (ix + (CG*smesh.nx+1) * iy)=0.0;
	  laplaceY        (ix + (CG*smesh.nx+1) * iy)=0.0;
	  momentum.GetVx()(lastnode - (ix + (CG*smesh.nx+1) * iy))=0.0;
	  momentum.GetVy()(lastnode - (ix + (CG*smesh.nx+1) * iy))=0.0;
	  laplaceX        (lastnode - (ix + (CG*smesh.nx+1) * iy))=0.0;
	  laplaceY        (lastnode - (ix + (CG*smesh.nx+1) * iy))=0.0;
	}
    Nextsim::VTK::write_cg_velocity(resultsdir + "/laplace", 0, momentum.GetVx(), momentum.GetVy(), smesh);
    Nextsim::VTK::write_cg_velocity(resultsdir + "/laplaceerror", 0, laplaceX, laplaceY, smesh);
    

    Nextsim::VTK::write_dg(resultsdir + "/e11", 0, momentum.GetE11(), smesh);
    Nextsim::VTK::write_dg(resultsdir + "/e11ex", 0, E11, smesh);

    
    std::cout << std::scientific << std::setprecision(3);
    std::cout << e11         << check(e11,test1[{CG,N}][0]) 
	      << "\t" << e12 << check(e12,test1[{CG,N}][1]) 
	      << "\t" << e22 << check(e22,test1[{CG,N}][2]) 
	      << "\t" << l1  << check(l1, test1[{CG,N}][3]) 
	      << "\t" << l2  << check(l2, test1[{CG,N}][4]) 
	      << std::endl;

}



template <int CG, int DGadvection>
void run_benchmark_laplace(const size_t N)
{
    //! Define the spatial mesh
  Nextsim::ParametricMesh smesh(Nextsim::SPHERICAL);
  create_mesh(smesh, 2*N, N);
  

    //! Compose name of output directory and create it
    std::string resultsdir = "TestGlobe_" + std::to_string(CG) + "_" + std::to_string(DGadvection) + "__" + std::to_string(N);
    std::filesystem::create_directory(resultsdir);
    
    //! Main class to handle the momentum equation. This class also stores the CG velocity vector
    Nextsim::CGParametricMomentum<CG> momentum(smesh, Nextsim::SPHERICAL);
 

    std::cout << CG << "\t" << N << "\t";
    
    ////////////////////////////////////////////////// Set Velocity
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetVx(), VX());
    Nextsim::Interpolations::Function2CG(smesh, momentum.GetVy(), VY());
    momentum.DirichletZero(momentum.GetVx());
    momentum.DirichletZero(momentum.GetVy());

    const Nextsim::CGVector<CG> vx_exact = momentum.GetVx();
    const Nextsim::CGVector<CG> vy_exact = momentum.GetVy();
    
    // disturb initial
    momentum.GetVx() *= 0.99;
    momentum.GetVy() *= 1.01;

    ////////////////////////////////////////////////// The rhs  f = -laplace
    Nextsim::CGVector<CG> laplaceX(smesh),laplaceY(smesh);
    Nextsim::Interpolations::Function2CG(smesh, laplaceX, LX());
    Nextsim::Interpolations::Function2CG(smesh, laplaceY, LY());
    Nextsim::VTK::write_cg_velocity(resultsdir + "/laplace", 0, laplaceX, laplaceY, smesh);

    // define the time mesh
    double ar = 4.0 * M_PI * SQR(Nextsim::EarthRadius)/(2.0*N*N);
    double dt = 0.05*ar/CG/CG;
    const double T = 5.0*Nextsim::EarthRadius*Nextsim::EarthRadius;
    const size_t NT = T/dt;

    // tmp vector for time-stepping
    Nextsim::CGVector<CG> tmpx(smesh), tmpy(smesh);
    
    for (size_t i=0; i<NT; ++i)
      {
	// compute strain
	momentum.ProjectCGVelocityToDGStrain();
	// laplace: Stress = Strain
	momentum.GetS11() = momentum.GetE11();
	momentum.GetS12() = momentum.GetE12();
	momentum.GetS22() = momentum.GetE22();

	// compute the divergence of the stress tensor
	momentum.DivergenceOfStress(1.0, tmpx, tmpy);
	
#pragma omp parallel for
	for (int i = 0; i < tmpx.rows(); ++i)
	  {
	    momentum.vx(i) += dt * (-laplaceX(i)+tmpx(i)/momentum.pmap.lumpedcgmass(i));
	    momentum.vy(i) += dt * (-laplaceY(i)+tmpy(i)/momentum.pmap.lumpedcgmass(i));
	  }
	momentum.DirichletZero(momentum.GetVx());
	momentum.DirichletZero(momentum.GetVy());
      }

    double l1 = sqrt(Nextsim::Interpolations::L2ErrorFunctionCG(smesh, momentum.vx, VX()))/SQR(Nextsim::EarthRadius);
    double l2 = sqrt(Nextsim::Interpolations::L2ErrorFunctionCG(smesh, momentum.vy, VY()))/SQR(Nextsim::EarthRadius);

    momentum.vx -= vx_exact;
    momentum.vy -= vy_exact;
    momentum.vx /= SQR(Nextsim::EarthRadius);
    momentum.vy /= SQR(Nextsim::EarthRadius);

    Nextsim::VTK::write_cg_velocity(resultsdir + "/velerror", 0, momentum.GetVx(), momentum.GetVy(), smesh);
    
    std::cout << std::scientific << std::setprecision(3);
    std::cout << l1 << check(l1,test2[{CG,N}][0]) << "\t"
	      << l2 << check(l2,test2[{CG,N}][1]) << std::endl;
}

int main()
{
  // ref values from 11/02/2023
  test1[{1,32}]  = {2.8347688048149036e-02,  4.6785274849483124e-02,  4.3911126871564121e-02,  3.3695047999915753e-02,  5.7941326833541695e-02};
  test1[{1,64}]  = {1.4158824215756512e-02,  2.3402756547335750e-02,  2.1930179939581396e-02,  8.5552566388180963e-03,  1.4698008059704821e-02};
  test1[{1,128}] = {7.0775287225742171e-03,  1.1702644278003502e-02,  1.0961898850099195e-02,  2.1519319208148672e-03,  3.6921445587101078e-03};
  test1[{2,32}]  = {3.5914636205682996e-04,  1.1807507404837177e-03,  1.2105698048397243e-03,  3.3835376857935368e-03,  4.3524556672044079e-03};
  test1[{2,64}]  = {8.9717798990046578e-05,  2.9532775918127228e-04,  3.0276932129230877e-04,  8.5461370907648313e-04,  1.0966014629455458e-03};
  test1[{2,128}] = {2.2425141391200837e-05,  7.3840681807791240e-05,  7.5700242049543274e-05,  2.1464766850264012e-04,  2.7488599304577381e-04};

  test2[{1,16}] = {5.2142793198142164e-03 ,4.0986719389934333e-03 };
  test2[{1,32}] = {1.2883627013626205e-03 ,1.0126458831291392e-03 };
  test2[{1,64}] = {3.2075639890935184e-04 ,2.5242056460174453e-04 };
  test2[{2,16}] = {1.2341037350435471e-04 ,1.6181065863764404e-04 };
  test2[{2,32}] = {1.5163120745444388e-05 ,1.9403500191718245e-05 };
  test2[{2,64}] = {2.0727071397833376e-06 ,2.3990803774422654e-06 };

  std::cout << "Convergence of the strain rate tensor and the Laplacian" << std::endl;
  std::cout << "CG\tN\tE11\t\tE12\t\tE22\t\tLAP1\t\tLAP2" << std::endl;
  for (auto NN : {32,64,128})
    run_benchmark_strain<1, 3>(NN);
  for (auto NN : {32,64,128})
    run_benchmark_strain<2, 3>(NN);


  std::cout << std::endl << "Convergence towards laplace solution" << std::endl;
  std::cout << "CG\tN\tLAP1\t\tLAP2" << std::endl;
  for (auto NN : {16, 32, 64})
    run_benchmark_laplace<1, 3>(NN);
  for (auto NN : {16, 32, 64})
    run_benchmark_laplace<2, 3>(NN);

}
