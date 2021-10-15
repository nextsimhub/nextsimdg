#include <cassert>
#include <iostream>
#include <chrono>

#include "dginitial.hpp"
#include "dglimiters.hpp"
#include "mesh.hpp"
#include "timemesh.hpp"
#include "dgvisu.hpp"
#include "dgtransport.hpp"

#define DG 2

int main()
{
  bool VTKOUTPUT = true;

  // Space mesh
  size_t N = 100;
  Nextsim::Mesh mesh(N, N, 1.0 / N);

  // main class to control the dg transport equation time stepping
  Nextsim::DGTransport<DG> dgtransport;

  // set the spacial mesh. dgtransport will adjust the size of the local vectors
  dgtransport.setmesh(mesh);

  // Initiallize the velocity field. Constant (in time) for this example
  Nextsim::CellVector<DG> vx(mesh), vy(mesh);

  Nextsim::L2ProjectInitial(mesh, vx, Nextsim::RotVXInitial());
  Nextsim::L2ProjectInitial(mesh, vy, Nextsim::RotVYInitial());

  // Vector to store the 'density' to be transported
  Nextsim::CellVector<DG> phi(mesh);
  Nextsim::L2ProjectInitial(mesh, phi, Nextsim::MixedInitial());

  Nextsim::Limiter<DG> limiter(mesh);

  // Time Mesh
  double TMAX = 0.1;
  double cfl = 0.1; // compute time step by cfl
  double k = cfl * mesh.h / 1.0;
  int NT = (static_cast<int>((TMAX / k + 1) / 100 + 1) * 100);

  Nextsim::TimeMesh timemesh(TMAX, NT);

  dgtransport.settimemesh(timemesh);

  if (VTKOUTPUT)
  {
    Nextsim::VTK::write_cellvector("Results/vx.vtk", vx, mesh);
    Nextsim::VTK::write_cellvector("Results/vy.vtk", vy, mesh);

    Nextsim::VTK::write_cellvector("Results/phi", 0, phi, mesh);
    Nextsim::VTK::write_dg("Results/dg", 0, phi, mesh);
  }

  int writestep = NT / 50; // for visualization

  //auto start = std::chrono::steady_clock::now();
  for (size_t iter = 1; iter <= timemesh.N; ++iter)
  {
    dgtransport.setvelocity(
        vx, vy); // sets the current velocity and averages it to edges
    dgtransport.step_heun(
        phi); // performs one time step with the 2nd Order Heun scheme

    if (iter % writestep == 0) // output & slope limiting
    {
      std::cout << "Time step " << iter << " / " << NT
                << "\t t = " << iter * timemesh.k << "\t" << phi.col(0).sum()
                << std::endl;

      if (VTKOUTPUT)
      {
        limiter.compute_vertex_based_limiter(phi);
        limiter.apply_vertex_based_limiter(phi);

        Nextsim::VTK::write_cellvector(
            "Results/phi", iter / writestep, phi, mesh);
        Nextsim::VTK::write_dg("Results/dg", iter / writestep, phi, mesh);
      }
    }
  }

  return 0;
}
