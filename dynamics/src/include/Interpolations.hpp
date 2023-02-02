/*!
 * @file Interpolations.hpp
 * @date Aug 9 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __INTERPOLATIONS_HPP
#define __INTERPOLATIONS_HPP

#include "cgVector.hpp"
#include "dgVector.hpp"

/*!
 * This file includes functions to interpolate between different vectors
 * that are used in the dynamics core, advection and momentum
 *
 * - Function2CG: Interpolates a function to a CG vector
 * - Function2DG: Projects a function to a DG vector
 * - CG2DG:       Projects a CG vector to a DG vector
 * - DG2GG:       Interpolates a DG Vector to a CG vector
 *
 * Interpolations just pick the node-wise values, where an interpolation from
 * DG to CG will average the values on the edges and the vertices
 *
 * Projections are performed as L2-projectio such that the local
 * DG mass matrix must in inverted.
 */

namespace Nextsim {

namespace Interpolations {

    //! interface class to provide an analytical function
    class Function {
    public:
        virtual double operator()(double x, double y) const = 0;
    };

    //! Interpolates an analytic function to a CG-Vector
    template <int CG>
    void Function2CG(const ParametricMesh& smesh, CGVector<CG>& dest, const Function& src);
    //! L2-Projection of an analytic function to a DG-Vector
    template <int DG>
    void Function2DG(const ParametricMesh& smesh, DGVector<DG>& dest, const Function& src, const COORDINATES CoordinateSystem);
    //! L2-Projection of CG-vector to a DG vector in Cartesian or Spherical coordinates
    template <int CG, int DG>
    void CG2DG(const ParametricMesh& smesh, DGVector<DG>& dest, const CGVector<CG>& src, const COORDINATES CoordinateSystem);
    //! Interpolation of DG-vector to a CG vector. Just averaging on edges / nodes
    template <int CG, int DG>
    void DG2CG(const ParametricMesh& smesh, CGVector<CG>& dest, const DGVector<DG>& src, const COORDINATES CoordinateSystem);

    //! Computes the L2 (integral) error between the DG-Vector and an analytic function
    template <int DG>
    double L2ErrorFunctionDG(const ParametricMesh& smesh, const DGVector<DG>& src, const Function& fct, const COORDINATES CoordinateSystem);

} /* namespace Interpolation */

} /* namespace Nextsim */

#endif /* __INTERPOLATIONS_HPP */
