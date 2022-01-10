/*----------------------------   cginitial.hpp     ---------------------------*/
#ifndef __cginitial_HPP
#define __cginitial_HPP
/*----------------------------   cginitial.hpp     ---------------------------*/

#include "cgvector.hpp"
#include "mesh.hpp"

namespace Nextsim {

//////////////////////////////////////////////////

//! Functions to project an analytical solution into the DG spaces

template <int CGdegree>
void InterpolateCG(const Mesh& mesh,
    CGVector<CGdegree>& phi,
    const InitialBase& initial);

template <>
void InterpolateCG(const Mesh& mesh,
    CGVector<2>& phi,
    const InitialBase& initial)
{
  assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == phi.rows());

#pragma omp parallel for
  for (size_t iy = 0; iy < 2*mesh.ny+1; ++iy) {
    const double Y = mesh.hy*0.5*iy;
    size_t ii = iy*(2*mesh.nx+1);
    
    for (size_t ix = 0; ix < 2*mesh.nx+1; ++ix, ++ii) {
      const double X = mesh.hx*0.5*ix;
      phi(ii) = initial(X,Y);
    }
  }
}

} // namespace Nextsim

/*----------------------------   cginitial.hpp     ---------------------------*/
/* end of #ifndef __cginitial_HPP */
#endif
/*----------------------------   cginitial.hpp     ---------------------------*/
