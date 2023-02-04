/*!
 * @file    NextsimDynamics.hpp
 * @date    Dec 13, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __NEXTSIMDYNAMICS_HPP
#define __NEXTSIMDYNAMICS_HPP


/*!
 * This file include all the definitions that are shared by nearly all nextsim-dynamics
 * files. 
 */

namespace Nextsim
{

  
  /*! 
   *
   * The different mappings used
   *
   * CARTESIAN: std. 2d cartesian map. No metric terms
   * SPHERICAL: mapping to spherical system in x=lon / y=lat
   *            ( R cos(lat)cos(lon), R cos(lat)sin(lon), R sin(lat) )
   */
  enum COORDINATES { SPHERICAL, CARTESIAN };

  /*!
   * Radius of the earth in [m]
   */
  constexpr double EarthRadius  = 6371000.0;  


  /*!
   * Computes the correct number of stress-unknowns
   * depending on the CG Discretization
   *
   * CG1 = 3, CG2 = 8
   */
#define CG2DGSTRESS(CG) ( (CG==1?3:(CG==2?8:-1) ) )

  
  // Number of local CG dofs per element
#define CGDOFS(CG) ( (CG==1)?4:9 )
  
  
  /*!
   * number of Gauss points depending on DG unknowns
   *  1 for dG(0)   DG=1
   * 2 for dG(1)   DG=3
   * 3 for dG(2)   DG=6 and DG=8
   */
#define GAUSSPOINTS(Q) (( (Q == 8) || (Q==6) ) ? 9 : (Q == 3) ? 4	\
			: (Q==1) ? 1 : -1)
#define GAUSSPOINTS1D(Q) (( (Q == 8) || (Q==6) ) ? 3 : (Q == 3) ? 2	\
			  : (Q==1) ? 1 : -1)

  /*!
   * The number of unknowns on an edge depending on DG, 
   * which is the number of unknowns in an element.
   */
#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))
  
  /*! 
   * the Float-Type used in NextSim. 
   * This is not yet fully checked or supported but might be required in 
   * for future use of GPU's
   */
  typedef double FloatType;

}


#endif // #define  __NEXTSIMDYNAMICS
