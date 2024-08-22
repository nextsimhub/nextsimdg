/*!
 * @file cgVector.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __CGVECTOR_HPP
#define __CGVECTOR_HPP

#include "ParametricMesh.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace Nextsim {

/*!
 * Stores a CG vector in std. Lagrangian Basis, e.g. having
 * (CGdegree*nx+1)*(CGdegree*ny+1) elements
 * Sorting lower left -> lower right -> ... -> upper right
 */
template <int CGdegree> class CGVector : public Eigen::Matrix<double, Eigen::Dynamic, 1> {
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> EigenCGVector;

    //! empty constructor
    CGVector() { }
    //! constructor setting size by mesh
    CGVector(const ParametricMesh& smesh)
        : EigenCGVector((CGdegree * smesh.nx + 1) * (CGdegree * smesh.ny + 1))
    {
    }

    //! resizes the vector and sets it to the mesh size
    void resize_by_mesh(const ParametricMesh& smesh)
    {
        EigenCGVector::resize((CGdegree * smesh.nx + 1) * (CGdegree * smesh.ny + 1));
    }

    // operations
    void zero() { EigenCGVector::setZero(); }

    // This method allows you to assign Eigen expressions to MyVectorType
    template <typename OtherDerived>
    CGVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, Eigen::Dynamic, 1>::operator=(other);
        return *this;
    }
};

} /* namespace Nextsim */

#endif /* __CGVECTOR_HPP */
