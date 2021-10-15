#include "dgtransport.hpp"
#include "dgtimestepping.hpp"

namespace Nextsim
{

  template<int DGdegree>
  void DGTransport<DGdegree>::setmesh(const Mesh &_mesh)
  {
    mesh = _mesh; // copy mesh
    
    tmp1.resize_by_mesh(mesh); // resize tmp-vectors for time stepping
    tmp2.resize_by_mesh(mesh);
    tmp3.resize_by_mesh(mesh);
    
    xvel_edgeY.resize_by_mesh(
			      mesh, EdgeType::Y); // resizes the velocity-edge-vectors
    yvel_edgeX.resize_by_mesh(mesh, EdgeType::X);
  }
  
  template<int DGdegree>
  void DGTransport<DGdegree>::settimemesh(const TimeMesh &_timemesh)
  {
    timemesh = _timemesh; // copy mesh
  }
  
  template<int DGdegree>
  void DGTransport<DGdegree>::setvelocity(const CellVector<DGdegree> &velx, const CellVector<DGdegree> &vely)
  {
    velxpointer = &velx; // copy pointers
    velypointer = &vely;
    
    // average the velocity to the edges
    average_to_edges_Y(mesh, xvel_edgeY, *velxpointer);
    average_to_edges_X(mesh, yvel_edgeX, *velypointer);
  }
  
  
  template<int DGdegree>
  void DGTransport<DGdegree>::step_fwdeuler(CellVector<DGdegree> &phi)
  {
    assert(velxpointer != NULL); // make sure that velocity is set
    assert(velypointer != NULL);
    
    fwdeuler_step<DGdegree>(mesh,
		      timemesh,
		      *velxpointer,
		      *velypointer,
		      xvel_edgeY,
		      yvel_edgeX,
		      phi,
		      tmp1);
    
    phi += tmp1;
  }

  template<int DGdegree>
  void DGTransport<DGdegree>::step_heun(CellVector<DGdegree> &phi)
  {
    assert(velxpointer != NULL); // make sure that velocity is set
    assert(velypointer != NULL);
    
    fwdeuler_step<DGdegree>(mesh,
		      timemesh,
		      *velxpointer,
		      *velypointer,
		      xvel_edgeY,
		      yvel_edgeX,
		      phi,
		      tmp1); // tmp1 = k * F(u)
    
    phi += tmp1; // phi = phi + k * F(u)     (i.e.: implicit Euler)
    
    fwdeuler_step<DGdegree>(mesh,
		      timemesh,
		      *velxpointer,
		      *velypointer,
		      xvel_edgeY,
		      yvel_edgeX,
		      phi,
		      tmp2); // tmp1 = k * F( u + k * F(u) )
    
    phi += 0.5 * (tmp2 - tmp1);
  }

  template class DGTransport<0>;
  template class DGTransport<1>;
  template class DGTransport<2>;
  

} // namespace Nextsim
