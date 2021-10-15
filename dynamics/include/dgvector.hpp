/*----------------------------   dgvector.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dgvector_H
#define __dgvector_H
/*----------------------------   dgvector.h     ---------------------------*/

#include <Eigen/Dense>
#include <iostream>

#include "mesh.hpp"

namespace Nextsim {
#define CELLDOFS(DGdegree) \
    (DGdegree == 0 ? 1 : (DGdegree == 1 ? 3 : (DGdegree == 2 ? 6 : -1)))

template <int DGdegree>
class LocalCellVector : public Eigen::Matrix<double, 1, CELLDOFS(DGdegree)> {
public:
    // required by Eigen
    LocalCellVector()
        : Eigen::Matrix<double, 1, CELLDOFS(DGdegree)>()
    {
    }

    template <typename OtherDerived>
    LocalCellVector(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix<double, 1, CELLDOFS(DGdegree)>(other)
    {
    }

    // This method allows you to assign Eigen expressions to MyVectorType
    template <typename OtherDerived>
    LocalCellVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, 1, CELLDOFS(DGdegree)>::operator=(other);
        return *this;
    }
};

template <int DGdegree>
class LocalEdgeVector : public Eigen::Matrix<double, 1, DGdegree + 1> {
public:
    LocalEdgeVector(void)
        : Eigen::Matrix<double, 1, DGdegree + 1>()
    {
    }

    template <typename OtherDerived>
    LocalEdgeVector(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix<double, 1, DGdegree + 1>(other)
    {
    }

    template <typename T1>
    LocalEdgeVector(const T1& t1)
        : Eigen::Matrix<double, 1, DGdegree + 1>(t1)
    {
    }
    template <typename T1, typename T2>
    LocalEdgeVector(const T1& t1, const T2& t2)
        : Eigen::Matrix<double, 1, DGdegree + 1>(t1, t2)
    {
    }
    template <typename T1, typename T2, typename T3>
    LocalEdgeVector(const T1& t1, const T2& t2, const T3& t3)
        : Eigen::Matrix<double, 1, DGdegree + 1>(t1, t2, t3)
    {
    }

    // This method allows you to assign Eigen expressions to MyVectorType
    template <typename OtherDerived>
    LocalEdgeVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->Eigen::Matrix<double, 1, DGdegree + 1>::operator=(other);
        return *this;
    }
};

/**
 * Stores coefficients of DGdegree vector
 *
 * Basis in local coordinate system:
 *
 * DGdegree 0:      1
 * DGdegree 1:      + (x-1/2), (y-1/2)
 * DGdegree 2:      + (x-1/2)^2-1/12, (y-1/2)^2-1/12, (x-1/2)(y-1/2)
 *
 **/
template <int DGdegree>
class CellVector
    : public Eigen::Matrix<double, Eigen::Dynamic, CELLDOFS(DGdegree)> {
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, CELLDOFS(DGdegree)>
        EigenCellVector;

    inline int dofs_in_cell() const { return CELLDOFS(DGdegree); }

    //! empty constructor
    CellVector() { }
    //! constructor setting size by mesh
    CellVector(const Mesh& mesh)
        : EigenCellVector(mesh.n, dofs_in_cell())
    {
    }

    //! resizes the vector and sets it to the mesh size
    void resize_by_mesh(const Mesh& mesh)
    {
        EigenCellVector::resize(mesh.n, dofs_in_cell());
    }

    // operations
    void zero() { EigenCellVector::setZero(); }
};

//! data set to store the type of the edges
typedef enum {
    none,
    X,
    Y
} EdgeType;

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
template <int DGdegree>
class EdgeVector : public Eigen::Matrix<double, Eigen::Dynamic, DGdegree + 1> {

public:
    //! Number of unknowns on each edge
    inline int dofs_in_edge() const { return DGdegree + 1; }

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
    EdgeVector(const Mesh& mesh, EdgeType et)
        : edgetype(et)
    {
        if (edgetype == X)
            Eigen::Matrix<double, Eigen::Dynamic, DGdegree + 1>::resize(
                mesh.nx * (mesh.ny + 1), dofs_in_edge());
        else if (edgetype == Y)
            Eigen::Matrix<double, Eigen::Dynamic, DGdegree + 1>::resize(
                (mesh.nx + 1) * mesh.ny, dofs_in_edge());
        else {
            std::cerr << "EdgeType must be set to X or Y" << std::endl;
            abort();
        }
    }

    //! resize vector
    void resize_by_mesh(const Mesh& mesh, EdgeType et)
    {
        edgetype = et;

        if (edgetype == X)
            Eigen::Matrix<double, Eigen::Dynamic, DGdegree + 1>::resize(
                mesh.nx * (mesh.ny + 1), dofs_in_edge());
        else if (edgetype == Y)
            Eigen::Matrix<double, Eigen::Dynamic, DGdegree + 1>::resize(
                (mesh.nx + 1) * mesh.ny, dofs_in_edge());
        else {
            std::cerr << "EdgeType must be set to X or Y" << std::endl;
            abort();
        }
    }

    // operations
    void zero()
    {
        Eigen::Matrix<double, Eigen::Dynamic, DGdegree + 1>::setZero();
    }
};

} // namespace Nextsim

/*----------------------------   dgvector.h     ---------------------------*/
/* end of #ifndef __dgvector_H */
#endif
/*----------------------------   dgvector.h     ---------------------------*/
