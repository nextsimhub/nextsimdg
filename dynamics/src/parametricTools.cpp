/*!
 * @file parametricTools.cpp
 * @date July 28, 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "ParametricTools.hpp"


namespace Nextsim {




    namespace ParametricTools
    {
  /*!
   * computes and fills the Q1/Q2 lumped mass matrix 
   */
  template<>
  void lumpedCGMassMatrix(const SasipMesh& smesh,
			  CGVector<1>& lumpedcgmass)
  {
    lumpedcgmass.resize_by_mesh(smesh);

#pragma omp parallel for
    for (size_t i=0;i<smesh.nnodes;++i)
      lumpedcgmass(i,0)=0;

    // parallelization not critical. called just once.

    const size_t sy = smesh.nx+1;
    size_t i=0;
    for (size_t iy=0;iy<smesh.ny;++iy)
      for (size_t ix=0;ix<smesh.nx;++ix,++i)
      {
	const double a = smesh.area(i);

	const int n0 = sy*iy+ix;
	lumpedcgmass(n0,0)+=0.25*a;
	lumpedcgmass(n0+1,0)+=0.25*a;
	lumpedcgmass(n0+sy,0)+=0.25*a;
	lumpedcgmass(n0+sy+1,0)+=0.25*a;
      }
  }

  template<>
  void lumpedCGMassMatrix(const SasipMesh& smesh,
			  CGVector<2>& lumpedcgmass)
  {
    lumpedcgmass.resize_by_mesh(smesh);

    for (size_t i=0;i<smesh.nnodes;++i)
      lumpedcgmass(i,0)=0;

    // parallelization not critical. called just once.
    
    const int sy = 2.0*smesh.nx+1;
    int i=0;
    for (size_t iy=0;iy<smesh.ny;++iy)
      for (size_t ix=0;ix<smesh.nx;++ix,++i)
	{
	const double a = smesh.area(i);
	const size_t n0 = 2*sy*iy+2*ix;
	lumpedcgmass(n0,0)+=a/36.;
	lumpedcgmass(n0+2,0)+=a/36.;
	lumpedcgmass(n0+2*sy,0)+=a/36.;
	lumpedcgmass(n0+2*sy+2,0)+=a/36.;
	
	lumpedcgmass(n0+1,0)+=a/9.;
	lumpedcgmass(n0+sy,0)+=a/9.;
	lumpedcgmass(n0+sy+2,0)+=a/9.;
	lumpedcgmass(n0+2*sy+1,0)+=a/9.;

	lumpedcgmass(n0+sy+1,0)+=a*4./9.;
	}

  }
    }


}
