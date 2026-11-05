/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICTOOLS_HPP
#define __PARAMETRICTOOLS_HPP

#include "NextsimDynamics.hpp"
#include "ParametricMesh.hpp"
#include "cgVector.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"
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
template <int CG, int DG> class SphericalTransformation {

public:
    /*!
     * These matrices realize the integration of (E, \grad phi) scaled with the
     * inverse mass matrix;
     */
    std::vector<Eigen::Matrix<FloatType, DG, CGDOFS(CG)>,
        Eigen::aligned_allocator<Eigen::Matrix<FloatType, DG, CGDOFS(CG)>>>
        iMgradX, iMgradY;

    std::vector<Eigen::Matrix<FloatType, CGDOFS(CG), DG>,
        Eigen::aligned_allocator<Eigen::Matrix<FloatType, CGDOFS(CG), DG>>>
        divS1, divS2;

    /*!
     * These matrices are M^-1 J w PSI_i(q)
     * Multiplied
     */
    std::vector<Eigen::Matrix<FloatType, DG, GAUSSPOINTS(DG)>,
        Eigen::aligned_allocator<Eigen::Matrix<FloatType, DG, GAUSSPOINTS(DG)>>>
        iMJwPSI;

    void BasicInit(const ParametricMesh& smesh);
};

namespace ParametricTools {

    /*!
     * computes and returns the gradient of the parametric map in the Gausspoints
     */
    template <int Q, typename T = FloatType>
    Eigen::Matrix<T, 2, Q * Q> dxT(const ParametricMesh& smesh, const size_t eid)
    {
        const Eigen::Matrix<T, 4, 2> coordinates = smesh.coordinatesOfElement<T>(eid);
        return coordinates.transpose() * PHIx<1, Q, T>;
    }
    template <int Q, typename T = FloatType>
    Eigen::Matrix<T, 2, Q * Q> dyT(const ParametricMesh& smesh, const size_t eid)
    {
        const Eigen::Matrix<T, 4, 2> coordinates = smesh.coordinatesOfElement<T>(eid);
        return coordinates.transpose() * PHIy<1, Q, T>;
    }

    /*!
     * computes and returns the degree of determinant of the transformation's Jacobian
     * depends on the number of gauss points Q
     */
    template <int Q, typename T = FloatType>
    Eigen::Matrix<T, 1, Q * Q> J(const ParametricMesh& smesh, const size_t eid)
    {
        // get the coordinates of the element as 4x2 - matrix
        const Eigen::Matrix<T, 4, 2> coordinates = smesh.coordinatesOfElement<T>(eid);
        const Eigen::Matrix<T, 2, Q * Q> dxT = coordinates.transpose() * PHIx<1, Q, T>;
        const Eigen::Matrix<T, 2, Q * Q> dyT = coordinates.transpose() * PHIy<1, Q, T>;

        // (dxT, dyT) is (dx T1, dx T2, dy T1, dy T2)
        return dxT.array().row(0) * dyT.array().row(1) - dxT.array().row(1) * dyT.array().row(0);
    }
    /*!
     * computes and returns the element mass matrix
     */
    template <int DG, typename T = FloatType>
    Eigen::Matrix<T, DG, DG> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
        constexpr int GP = GAUSSPOINTS1D(DG);
        return (PSI<DG, GP, T>.array().rowwise()
                   * (GAUSSWEIGHTS<GP, T>.array() * J<GP, T>(smesh, eid).array()))
                   .matrix()
            * PSI<DG, GP, T>.transpose();
    }

    /*!
     * computes and retunrs the coordinates of the NGP^2 gauss points
     * in the physical element with index eid
     */
    template <int NGP1d, typename T = FloatType>
    Eigen::Matrix<T, 2, NGP1d * NGP1d> getGaussPointsInElement(
        const ParametricMesh& smesh, const size_t eid)
    {
        return smesh.coordinatesOfElement<T>(eid).transpose() * PHI<1, NGP1d, T>;
    }

    /*!
     * computes and returns the coordinates of the NGP gauss points
     * in the physical element along the edge with index eid
     */
    template <int NGP1d>
    inline Eigen::Matrix<FloatType, 2, NGP1d> getGaussPointsOnEdgeX(
        const ParametricMesh& smesh, const size_t eid)
    {
        return smesh.coordinatesOfEdgeX(eid).transpose() * PHI1d<1, NGP1d>;
    }
    template <int NGP1d>
    inline Eigen::Matrix<FloatType, 2, NGP1d> getGaussPointsOnEdgeY(
        const ParametricMesh& smesh, const size_t eid)
    {
        return smesh.coordinatesOfEdgeY(eid).transpose() * PHI1d<1, NGP1d>;
    }

}

namespace SphericalTools {
    /*!
     * computes and returns the element mass matrix
     */
    template <int DG, typename T = FloatType>
    inline Eigen::Matrix<T, DG, DG> massMatrix(const ParametricMesh& smesh, const size_t eid)
    {
        return (PSI<DG, GAUSSPOINTS1D(DG), T>.array().rowwise()
                   * (GAUSSWEIGHTS<GAUSSPOINTS1D(DG), T>.array()
                       * ParametricTools::J<GAUSSPOINTS1D(DG), T>(smesh, eid).array()
                       * (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG), T>(smesh, eid)
                               .row(1)
                               .array())
                             .cos()))
                   .matrix()
            * PSI<DG, GAUSSPOINTS1D(DG), T>.transpose();
    }

}

} /* namespace Nextsim */

#endif /* __PARAMETRICTOOLS_HPP */
