/*!
 * @file    ParametricMap.hpp
 * @date    Dec 13, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICMAP_HPP
#define __PARAMETRICMAP_HPP

#include "ParametricMesh.hpp"
#include "NextsimDynamics.hpp"

namespace Nextsim
{
  
  template<int DG>
  class ParametricMap
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
    

    ParametricMap(const ParametricMesh& sm, COORDINATES ty) : smesh(sm), type(ty)
    {}

    //! initialization of the different forms

    //! For cell_term in Advection. Fills AdvectionCellTermX/Y
    void InitializeAdvectionCellTerms();
    //! For cell_term in Advection. Fills AdvectionCellTermX/Y
    void InitializeInverseDGMassMatrix();
     
  };
}


#endif // #define  __PARAMETRICMAP_HPP
