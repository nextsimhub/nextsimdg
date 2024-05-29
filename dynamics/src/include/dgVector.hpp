/*!
 * @file dgVector.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGVECTOR_HPP
#define __DGVECTOR_HPP
#include "ParametricMesh.hpp"

#include <Eigen/Dense>
#include <iostream>

namespace Nextsim {
// #define CELLDOFS(DGdegree) (DGdegree == 0 ? 1 : (DGdegree == 1 ? 3 : (DGdegree == 2 ? 6 : -1)))

template <int DG> class LocalDGVector : public Eigen::Matrix<double, 1, DG> {
public:
    // required by Eigen
    LocalDGVector()
        : Eigen::Matrix<double, 1, DG>()
    {
    }

    template <typename OtherDerived>
    LocalDGVector(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix<double, 1, DG>(other)
    {
    }

    // This method allows you to assign Eigen expressions to MyVectorType
    template <typename OtherDerived>
    LocalDGVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, 1, DG>::operator=(other);
        return *this;
    }
};

template <int DG> class LocalEdgeVector : public Eigen::Matrix<double, 1, DG> {
public:
    LocalEdgeVector(void)
        : Eigen::Matrix<double, 1, DG>()
    {
    }

    template <typename OtherDerived>
    LocalEdgeVector(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix<double, 1, DG>(other)
    {
    }

    template <typename T1>
    LocalEdgeVector(const T1& t1)
        : Eigen::Matrix<double, 1, DG>(t1)
    {
    }
    template <typename T1, typename T2>
    LocalEdgeVector(const T1& t1, const T2& t2)
        : Eigen::Matrix<double, 1, DG>(t1, t2)
    {
    }
    template <typename T1, typename T2, typename T3>
    LocalEdgeVector(const T1& t1, const T2& t2, const T3& t3)
        : Eigen::Matrix<double, 1, DG>(t1, t2, t3)
    {
    }

    // This method allows you to assign Eigen expressions to MyVectorType
    template <typename OtherDerived>
    LocalEdgeVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, 1, DG>::operator=(other);
        return *this;
    }
};

/**
 * Stores coefficients of DGdegree vector
 *
 * Basis in local coordinate system:
 *
 * DGdegree 0:      1
 * DGdegree 1-2:    + (x-1/2),
 * DGdegree 3-5:    + (x-1/2)^2-1/12, (y-1/2)^2-1/12, (x-1/2)(y-1/2)
 * DGdegree 6-7:    + (y-1/2)(x-1/2)^2-1/12, (x-1/2)(y-1/2)^2-1/12
 *
 **/
template <int DG>
class DGVector : public Eigen::Matrix<double, Eigen::Dynamic, DG,
                     (DG == 1) ? Eigen::ColMajor : Eigen::RowMajor> {
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, DG, (DG == 1) ? Eigen::ColMajor : Eigen::RowMajor>
        EigenDGVector;

    inline int dofs_in_cell() const { return DG; }

    //! empty constructor
    DGVector() { }
    //! constructor setting size by mesh
    DGVector(const ParametricMesh& smesh)
        : EigenDGVector(smesh.nelements, dofs_in_cell())
    {
    }

    //! resizes the vector and sets it to the mesh size
    void resize_by_mesh(const ParametricMesh& smesh)
    {
        EigenDGVector::resize(smesh.nelements, dofs_in_cell());
    }

    // operations
    void zero() { EigenDGVector::setZero(); }

    // This method allows you to assign Eigen expressions to MyVectorType
    template <typename OtherDerived>
    DGVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, Eigen::Dynamic, DG,
            (DG == 1) ? Eigen::ColMajor : Eigen::RowMajor>::operator=(other);
        return *this;
    }
    template <typename OtherDerived>
    DGVector& operator+=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, Eigen::Dynamic, DG,
            (DG == 1) ? Eigen::ColMajor : Eigen::RowMajor>::operator+=(other);
        return *this;
    }
};

//! data set to store the type of the edges
typedef enum { none, X, Y } EdgeType;

/*!
 * Stores coefficients of DGdegree vector on edges
 *
 * DGdegree 0:    1
 * DGdegree 1: +  (t-1/2)
 * DGdegree 2: +  (t-1/2)^2-1/12
 *
 * edges (t) are oriented in positive x- and y-direction
 *
 *
 **/
template <int DG> class EdgeVector : public Eigen::Matrix<double, Eigen::Dynamic, DG> {

public:
    //! Number of unknowns on each edge
    inline int dofs_in_edge() const { return DG; }

    /*!
     * Type of the edge:
     * none: not assigned
     * X:    edge is parallel to X-axes. Vector has size nx * (ny+1)
     * Y:    edge is parallel to Y-axes. Vector has size (nx+1) * ny
     */

    EdgeType edgetype;

    // empty constructor
    EdgeVector()
        : edgetype(none)
    {
    }
    //! constructor setting size by mesh
    EdgeVector(const ParametricMesh& smesh, EdgeType et)
        : edgetype(et)
    {
        if (edgetype == X)
            Eigen::Matrix<double, Eigen::Dynamic, DG>::resize(
                smesh.nx * (smesh.ny + 1), dofs_in_edge());
        else if (edgetype == Y)
            Eigen::Matrix<double, Eigen::Dynamic, DG>::resize(
                (smesh.nx + 1) * smesh.ny, dofs_in_edge());
        else {
            std::cerr << "EdgeType must be set to X or Y" << std::endl;
            abort();
        }
    }

    //! resize vector
    void resize_by_mesh(const ParametricMesh& smesh, EdgeType et)
    {
        edgetype = et;

        if (edgetype == X)
            Eigen::Matrix<double, Eigen::Dynamic, DG>::resize(
                smesh.nx * (smesh.ny + 1), dofs_in_edge());
        else if (edgetype == Y)
            Eigen::Matrix<double, Eigen::Dynamic, DG>::resize(
                (smesh.nx + 1) * smesh.ny, dofs_in_edge());
        else {
            std::cerr << "EdgeType must be set to X or Y" << std::endl;
            abort();
        }
    }

    // operations
    void zero() { Eigen::Matrix<double, Eigen::Dynamic, DG>::setZero(); }
};

} /* namespace Nextsim */

#endif /* __DGVECTOR_HPP */
