/*!
 * @file    ParametricTools.hpp
 * @date    July 10, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICTOOLS_HPP
#define __PARAMETRICTOOLS_HPP

#include "SasipMesh.hpp"
#include "codeGenerationDGinGauss.hpp"

/*!
 * Selection of functions required to do integration on the parametric Sasip-Mesh
 *
 * - massMatrix<DG>(eid)     returns the element mass matrix on mesh elemnt id with DG dofs
 *
 */

namespace Nextsim {

  namespace ParametricTools
  {

    /*!
     * computes and returns the gradient of the parametric map in the Gausspoints
     */
    template<int Q>
    Eigen::Matrix<Nextsim::FloatType, 2, Q*Q> dxT(const SasipMesh& smesh, const size_t eid);
    template<>
    Eigen::Matrix<Nextsim::FloatType, 2, 1> dxT<1>(const SasipMesh& smesh, const size_t eid)
    {
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);
      return coordinates.transpose() * CG_Q1x_1;
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 2, 4> dxT<2>(const SasipMesh& smesh, const size_t eid)
    {
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);
      return coordinates.transpose() * CG_Q1x_2;
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 2, 9> dxT<3>(const SasipMesh& smesh, const size_t eid)
    {
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);
      return coordinates.transpose() * CG_Q1x_3;
    }
    template<int Q>
    Eigen::Matrix<Nextsim::FloatType, 2, Q*Q> dyT(const SasipMesh& smesh, const size_t eid);
    template<>
    Eigen::Matrix<Nextsim::FloatType, 2, 1> dyT<1>(const SasipMesh& smesh, const size_t eid)
    {
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);
      return coordinates.transpose() * CG_Q1y_1;
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 2, 4> dyT<2>(const SasipMesh& smesh, const size_t eid)
    {
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);
      return coordinates.transpose() * CG_Q1y_2;
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 2, 9> dyT<3>(const SasipMesh& smesh, const size_t eid)
    {
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);
      return coordinates.transpose() * CG_Q1y_3;
    }


    /*!
     * computes and returns the degree of determinant of the transformation's Jacobian
     * depends on the number of gauss points Q
     */
    template<int Q>
    Eigen::Matrix<Nextsim::FloatType, 1, Q*Q> J(const SasipMesh& smesh, const size_t eid);

    template<>
    Eigen::Matrix<Nextsim::FloatType, 1, 1> J<1>(const SasipMesh& smesh, const size_t eid)
    {
      // get the coordinates of the element as 4x2 - matrix
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);

      const Eigen::Matrix<Nextsim::FloatType, 2,1> dxT = coordinates.transpose() * CG_Q1x_1;
      const Eigen::Matrix<Nextsim::FloatType, 2,1> dyT = coordinates.transpose() * CG_Q1y_1;

      // (dxT, dyT) is (dx T1, dx T2, dy T1, dy T2)
      return dxT.array().row(0)*dyT.array().row(1)-dxT.array().row(1)*dyT.array().row(0);
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 1, 4> J<2>(const SasipMesh& smesh, const size_t eid)
    {
      // get the coordinates of the element as 4x2 - matrix
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);

      const Eigen::Matrix<Nextsim::FloatType, 2,4> dxT = coordinates.transpose() * CG_Q1x_2;
      const Eigen::Matrix<Nextsim::FloatType, 2,4> dyT = coordinates.transpose() * CG_Q1y_2;

      // (dxT, dyT) is (dx T1, dx T2, dy T1, dy T2)
      return dxT.array().row(0)*dyT.array().row(1)-dxT.array().row(1)*dyT.array().row(0);
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 1, 9> J<3>(const SasipMesh& smesh, const size_t eid)
    {
      // get the coordinates of the element as 4x2 - matrix
      const Eigen::Matrix<Nextsim::FloatType, 4,2> coordinates
	= smesh.coordinatesOfElement(eid);

      const Eigen::Matrix<Nextsim::FloatType, 2,9> dxT = coordinates.transpose() * CG_Q1x_3;
      const Eigen::Matrix<Nextsim::FloatType, 2,9> dyT = coordinates.transpose() * CG_Q1y_3;

      // (dxT, dyT) is (dx T1, dx T2, dy T1, dy T2)
      return dxT.array().row(0)*dyT.array().row(1)-dxT.array().row(1)*dyT.array().row(0);
    }

    /*! 
     * computes and returns the element mass matrix
     */ 
    template<int DG>
    Eigen::Matrix<Nextsim::FloatType, DG, DG> massMatrix(const SasipMesh& smesh, const size_t eid);

    template<>
    Eigen::Matrix<Nextsim::FloatType, 1, 1> massMatrix(const SasipMesh& smesh, const size_t eid)
    {
      return Eigen::Matrix<Nextsim::FloatType, 1,1> (smesh.area(eid));

      // mit 1 GP. Reicht das???
      // return (BiG31.array().rowwise() * (GAUSSWEIGHTS_1.array() * J<1>(smesh,eid).array())).matrix() * BiG31.transpose();
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 3, 3> massMatrix(const SasipMesh& smesh, const size_t eid)
    {
      return (BiG32.array().rowwise() * (GAUSSWEIGHTS_2.array() * J<2>(smesh,eid).array())).matrix() * BiG32.transpose();

      // mit 1 GP. Reicht das???
      // return (BiG31.array().rowwise() * (GAUSSWEIGHTS_1.array() * J<1>(smesh,eid).array())).matrix() * BiG31.transpose();
    }
    template<>
    Eigen::Matrix<Nextsim::FloatType, 6, 6> massMatrix(const SasipMesh& smesh, const size_t eid)
    {
      return (BiG63.array().rowwise() * (GAUSSWEIGHTS_3.array() * J<3>(smesh,eid).array())).matrix() * BiG63.transpose();
    }


    /*!
     * computes and returs the coordinates of the Q^2 gauss points
     * in the physical element with index eid
     */
    Eigen::Matrix<Nextsim::FloatType, 2, 1> getGaussPointsInElement1(const SasipMesh& smesh, const size_t eid)
    {
      return smesh.coordinatesOfElement(eid).transpose() * CG_Q1_1;
    }
    Eigen::Matrix<Nextsim::FloatType, 2, 4> getGaussPointsInElement2(const SasipMesh& smesh, const size_t eid)
    {
      return smesh.coordinatesOfElement(eid).transpose() * CG_Q1_2;
    }
    Eigen::Matrix<Nextsim::FloatType, 2, 9> getGaussPointsInElement3(const SasipMesh& smesh, const size_t eid)
    {
      return smesh.coordinatesOfElement(eid).transpose() * CG_Q1_3;
    }
    

    
  }

} /* namespace Nextsim */

#endif /* __PARAMETRICTOOLS_HPP */
