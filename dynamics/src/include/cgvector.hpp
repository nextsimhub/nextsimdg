/*----------------------------   cgvector.hpp     ---------------------------*/
#ifndef __cgvector_HPP
#define __cgvector_HPP
/*----------------------------   cgvector.hpp     ---------------------------*/

/*!
 * @file   cgvector.hpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include <Eigen/Dense>
#include <iostream>

#include "mesh.hpp"

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
    CGVector(const Mesh& mesh)
        : EigenCGVector((CGdegree * mesh.nx + 1) * (CGdegree * mesh.ny + 1))
    {
    }

    //! resizes the vector and sets it to the mesh size
    void resize_by_mesh(const Mesh& mesh)
    {
        EigenCGVector::resize((CGdegree * mesh.nx + 1) * (CGdegree * mesh.ny + 1));
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

}

/*----------------------------   cgvector.hpp     ---------------------------*/
/* end of #ifndef __cgvector_HPP */
#endif
/*----------------------------   cgvector.hpp     ---------------------------*/
