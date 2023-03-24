/*!
 * @file    ParametricTools.hpp
 * @date    July 10, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICTOOLS_HPP
#define __PARAMETRICTOOLS_HPP

#include "NextsimDynamics.hpp"
#include "ParametricMesh.hpp"
#include "cgVector.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "codeGenerationCGinGauss.hpp"
#include <Eigen/StdVector>
#include <vector>

/*!
 * Selection of functions required to do integration on the parametric Sasip-Mesh
 *
 * - massMatrix<DG>(eid)     returns the element mass matrix on mesh elemnt id with DG dofs
 *
 */

namespace Nextsim {


/*!
 * This class stores precomputed values (matrices in each mesh element)
 * that are required to efficiently perform the advection scheme and the
 * momentum equation on the parametric mesh.
 *
 * The values must be initialized once for a mesh. The storage requirements
 * are substantial. An alternative would be to recompute the quantities
 * whenever required. (see ParametricTools)
 *
 * Additional storage per element for CG2 / DG8
 * 1) divS1 + divS2     =  2 * 9*8 = 144
 * 2) iMgradX + iMgradY =  2 * 9*8 = 144
 * 3) iMJwPSI           =      8*9 = 72
 * 360 doubles = 2.88 kB
 */
template <int CG, int DG>
class SphericalTransformation {

public:

    /*!
     * These matrices realize the integration of (E, \grad phi) scaled with the
     * inverse mass matrix;
     */
    std::vector<Eigen::Matrix<Nextsim::FloatType, DG, CGDOFS(CG)>,
        Eigen::aligned_allocator<Eigen::Matrix<Nextsim::FloatType, DG, CGDOFS(CG)>>>
        iMgradX, iMgradY;

    std::vector<Eigen::Matrix<Nextsim::FloatType, CGDOFS(CG), DG>,
        Eigen::aligned_allocator<Eigen::Matrix<Nextsim::FloatType, CGDOFS(CG), DG>>>
        divS1, divS2;

    /*!
     * These matrices are M^-1 J w PSI_i(q)
     * Multiplied
     */
    std::vector<Eigen::Matrix<Nextsim::FloatType, DG, GAUSSPOINTS(DG)>,
        Eigen::aligned_allocator<Eigen::Matrix<Nextsim::FloatType, DG, GAUSSPOINTS(DG)>>>
        iMJwPSI;

    void BasicInit(const ParametricMesh& smesh);
};
  
namespace ParametricTools {

    /*!
     * computes and returns the gradient of the parametric map in the Gausspoints
     */
    template <int Q>
    inline Eigen::Matrix<Nextsim::FloatType, 2, Q * Q> dxT(const ParametricMesh& smesh, const size_t eid)
    {
        const Eigen::Matrix<Nextsim::FloatType, 4, 2> coordinates
            = smesh.coordinatesOfElement(eid);
        return coordinates.transpose() * PHIx<1,Q>;
    }
    template <int Q>
    inline Eigen::Matrix<Nextsim::FloatType, 2, Q * Q> dyT(const ParametricMesh& smesh, const size_t eid)
    {
        const Eigen::Matrix<Nextsim::FloatType, 4, 2> coordinates
            = smesh.coordinatesOfElement(eid);
        return coordinates.transpose() * PHIy<1,Q>;
    }

    /*!
     * computes and returns the degree of determinant of the transformation's Jacobian
     * depends on the number of gauss points Q
     */
    template <int Q>
    inline Eigen::Matrix<Nextsim::FloatType, 1, Q * Q> J(const ParametricMesh& smesh, const size_t eid)
    {
      // get the coordinates of the element as 4x2 - matrix
        const Eigen::Matrix<Nextsim::FloatType, 4, 2> coordinates
            = smesh.coordinatesOfElement(eid);
        const Eigen::Matrix<Nextsim::FloatType, 2, Q*Q> dxT = coordinates.transpose() * PHIx<1,Q>;
        const Eigen::Matrix<Nextsim::FloatType, 2, Q*Q> dyT = coordinates.transpose() * PHIy<1,Q>;

        // (dxT, dyT) is (dx T1, dx T2, dy T1, dy T2)
        return dxT.array().row(0) * dyT.array().row(1) - dxT.array().row(1) * dyT.array().row(0);
    }
    /*!
     * computes and returns the element mass matrix
     */
    template <int DG>
    inline Eigen::Matrix<Nextsim::FloatType, DG, DG> massMatrix(const ParametricMesh& smesh, const size_t eid);

    template <>
    inline Eigen::Matrix<Nextsim::FloatType, 1, 1> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
        return Eigen::Matrix<Nextsim::FloatType, 1, 1>(smesh.area(eid));

        // mit 1 GP. Reicht das???
        // return (PSI31.array().rowwise() * (GAUSSWEIGHTS_1.array() * J<1>(smesh,eid).array())).matrix() * PSI31.transpose();
    }
    template <>
    inline Eigen::Matrix<Nextsim::FloatType, 3, 3> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
        return (PSI<3, 2>.array().rowwise() * (GAUSSWEIGHTS<2>.array() * J<2>(smesh, eid).array())).matrix() * PSI<3, 2>.transpose();

        // mit 1 GP. Reicht das???
        // return (PSI31.array().rowwise() * (GAUSSWEIGHTS_1.array() * J<1>(smesh,eid).array())).matrix() * PSI31.transpose();
    }
    template <>
    inline Eigen::Matrix<Nextsim::FloatType, 6, 6> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
        return (PSI<6, 3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J<3>(smesh, eid).array())).matrix() * PSI<6, 3>.transpose();
    }
    template <>
    inline Eigen::Matrix<Nextsim::FloatType, 8, 8> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
        return (PSI<8, 3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J<3>(smesh, eid).array())).matrix() * PSI<8, 3>.transpose();
    }

    /*!
     * computes and retunrs the coordinates of the NGP^2 gauss points
     * in the physical element with index eid
     */
  template<int NGP1d>
    inline Eigen::Matrix<Nextsim::FloatType, 2, NGP1d*NGP1d> getGaussPointsInElement(const ParametricMesh& smesh, const size_t eid)
    {
      return smesh.coordinatesOfElement(eid).transpose() * PHI<1,NGP1d>;
    }

  /*!
   * computes and returns the coordinates of the NGP gauss points
     * in the physical element along the edge with index eid
     */
  template<int NGP1d>
    inline Eigen::Matrix<Nextsim::FloatType, 2, NGP1d> getGaussPointsOnEdgeX(const ParametricMesh& smesh, const size_t eid)
    {
      return smesh.coordinatesOfEdgeX(eid).transpose() * PHI1d<1,NGP1d>;
    }
  template<int NGP1d>
    inline Eigen::Matrix<Nextsim::FloatType, 2, NGP1d> getGaussPointsOnEdgeY(const ParametricMesh& smesh, const size_t eid)
    {
      return smesh.coordinatesOfEdgeY(eid).transpose() * PHI1d<1,NGP1d>;
    }


}


namespace SphericalTools {
    /*!
     * computes and returns the element mass matrix
     */
    template <int DG>
    inline Eigen::Matrix<Nextsim::FloatType, DG, DG> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
      return (PSI<DG, GAUSSPOINTS1D(DG) >.array().rowwise() *
	      (GAUSSWEIGHTS<GAUSSPOINTS1D(DG) >.array() *
	       ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, eid).array() *
	       (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(smesh, eid).row(1).array()).cos()
	       )).matrix() * PSI<DG, GAUSSPOINTS1D(DG)>.transpose();
    }

}
  

} /* namespace Nextsim */

#endif /* __PARAMETRICTOOLS_HPP */
