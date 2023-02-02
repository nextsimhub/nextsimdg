/*!
 * @file    ParametricMap.hpp
 * @date    Dec 13, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICMAP_HPP
#define __PARAMETRICMAP_HPP

#include "ParametricMesh.hpp"
#include "NextsimDynamics.hpp"
#include "cgVector.hpp"

namespace Nextsim
{

  /*!
   * Stores precomputed matrices and stencils that are required
   * for the advection. 
   *
   * The coordiante system is encoded into the matrices
   */
  template<int DG>
  class ParametricTransportMap
  {
    //! Reference to the map. Given with constructor
    const ParametricMesh& smesh;

    const COORDINATES type;

    //! What type of 
  public:

    //! These terms are required for the cell-term in the advection  -(vA, nabla PHI)
    std::vector< Eigen::Matrix<Nextsim::FloatType, DG, GP(DG)*GP(DG) > > AdvectionCellTermX;
    std::vector< Eigen::Matrix<Nextsim::FloatType, DG, GP(DG)*GP(DG) > > AdvectionCellTermY;

    //! The inverse of the dG mass matrix
    std::vector< Eigen::Matrix<Nextsim::FloatType, DG, DG> > InverseDGMassMatrix;
    

    ParametricTransportMap(const ParametricMesh& sm, COORDINATES ty) : smesh(sm), type(ty)
    {}

    //! initialization of the different forms

    //! For cell_term in Advection. Fills AdvectionCellTermX/Y
    void InitializeAdvectionCellTerms();

    
    /*!
     *
     *  For cell_term in Advection. Fills AdvectionCellTermX/Y
     *
     * In Spherical coordinates the inverse mass is scaled by 1/Radius to account for
     * the proper scaling. 
     *
     * R^2 (v',phi) - R (vA, nabla phi) + R <va * N, phi> = 0
     */
    void InitializeInverseDGMassMatrix();
     
  };




  /*!
   * Stores precomputed matrices and stencils that are required
   * for the momentum equation
   *
   * The coordiante system is encoded into the matrices
   */
  template<int CG>
  class ParametricMomentumMap
  {
    //! Reference to the map. Given with constructor
    const ParametricMesh& smesh;

    //! What type of coordinate system is used
    const COORDINATES type;

  public:

    //! Vector to store the lumpes mass matrix. Is directly initialized when the mesh is known
    CGVector<CG> lumpedcgmass;

    /*!
     * These matrices realize the integration of (-div S, phi) = (S, nabla phi)
     * as matrix-vector producs (divS1 * S11 + divS2 * S12 ; divS1 * S21 + divS2 * S22)
     * [ where S12= S21 ]
     */
    std::vector<Eigen::Matrix<Nextsim::FloatType, CGDOFS(CG),CG2DGSTRESS(CG)>,
		Eigen::aligned_allocator<Eigen::Matrix<Nextsim::FloatType, CGDOFS(CG), CG2DGSTRESS(CG)>>>
    divS1, divS2;
    
    
    ParametricMomentumMap(const ParametricMesh& sm, COORDINATES ty) : smesh(sm), type(ty)
    {}

    //! initialization of the different forms


    //! initializes the lumped mass matrix for the stress update
    void InitializeLumpedCGMassMatrix();
    //! initializes div-matrices for the stress update
    void InitializeDivSMatrices();
  };



}


#endif // #define  __PARAMETRICMAP_HPP
