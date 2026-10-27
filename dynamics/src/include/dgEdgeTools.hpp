/*!
 * @file dgEdgeTools.hpp
 *
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 *
 * Restriction of the DG coefficients of a cell to the coefficients of the
 * trace on one of its edges, expressed in the 1D edge basis (EDGEDOFS(DG)
 * coefficients). Used by the DG transport and the LDG momentum.
 */

#ifndef __DGEDGETOOLS_HPP
#define __DGEDGETOOLS_HPP

#include "NextsimDynamics.hpp"
#include "dgVector.hpp"

namespace Nextsim {

//! returns the localization of the cell vector to the edges
/*!
 * writes the cell-basis on the edges in the edge basis
 *
 * CELL: (X = x-1/2, Y = y-1/2)
 * index 0:      1
 * index 1-2:    X, Y
 * index 3:      XY
 * index 4-5:    X^2-1/12, Y^2-1/12
 * index 6-7:    Y(X^2-1/12), X(Y^2-1/12)
 * index 8:      (X^2-1/12)(Y^2-1/12)
 *
 * EDGE:
 * DGdegree 0:    1
 * DGdegree 1: +  (t-1/2)
 * DGdegree 2: +  (t-1/2)^2-1/12
 *
 */
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> leftedgeofcell(const DGVector<DG>& cv, size_t eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> rightedgeofcell(const DGVector<DG>& cv, size_t eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> bottomedgeofcell(const DGVector<DG>& cv, size_t eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> topedgeofcell(const DGVector<DG>& cv, size_t eid);

// dG0 (1 in cell, 1 on edge)
template <> inline Eigen::Matrix<FloatType, 1, 1> leftedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}
template <> inline Eigen::Matrix<FloatType, 1, 1> rightedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}
template <>
inline Eigen::Matrix<FloatType, 1, 1> bottomedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}
template <> inline Eigen::Matrix<FloatType, 1, 1> topedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}

// dG1 (3 in cell, 2 on edge)
template <> inline Eigen::Matrix<FloatType, 1, 2> leftedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5_ft * cv(eid, 1), cv(eid, 2));
}
template <> inline Eigen::Matrix<FloatType, 1, 2> rightedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5_ft * cv(eid, 1), cv(eid, 2));
}
template <>
inline Eigen::Matrix<FloatType, 1, 2> bottomedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5_ft * cv(eid, 2), cv(eid, 1));
}
template <> inline Eigen::Matrix<FloatType, 1, 2> topedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5_ft * cv(eid, 2), cv(eid, 1));
}

// dQ1 (4 in cell, 2 on edge)
template <> inline Eigen::Matrix<FloatType, 1, 2> leftedgeofcell(const DGVector<4>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(
        cv(eid, 0) - 0.5_ft * cv(eid, 1), cv(eid, 2) - 0.5_ft * cv(eid, 3));
}
template <> inline Eigen::Matrix<FloatType, 1, 2> rightedgeofcell(const DGVector<4>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(
        cv(eid, 0) + 0.5_ft * cv(eid, 1), cv(eid, 2) + 0.5_ft * cv(eid, 3));
}
template <>
inline Eigen::Matrix<FloatType, 1, 2> bottomedgeofcell(const DGVector<4>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(
        cv(eid, 0) - 0.5_ft * cv(eid, 2), cv(eid, 1) - 0.5_ft * cv(eid, 3));
}
template <> inline Eigen::Matrix<FloatType, 1, 2> topedgeofcell(const DGVector<4>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(
        cv(eid, 0) + 0.5_ft * cv(eid, 2), cv(eid, 1) + 0.5_ft * cv(eid, 3));
}

// dG2 (6 in cell, 3 on edge)
template <> inline Eigen::Matrix<FloatType, 1, 3> leftedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5_ft * cv(eid, 1) + 1. / 6. * cv(eid, 4),
        cv(eid, 2) - 0.5_ft * cv(eid, 3), cv(eid, 5));
}
template <> inline Eigen::Matrix<FloatType, 1, 3> rightedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5_ft * cv(eid, 1) + 1. / 6. * cv(eid, 4),
        cv(eid, 2) + 0.5_ft * cv(eid, 3), cv(eid, 5));
}
template <>
inline Eigen::Matrix<FloatType, 1, 3> bottomedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5_ft * cv(eid, 2) + 1. / 6. * cv(eid, 5),
        cv(eid, 1) - 0.5_ft * cv(eid, 3), cv(eid, 4));
}
template <> inline Eigen::Matrix<FloatType, 1, 3> topedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5_ft * cv(eid, 2) + 1. / 6. * cv(eid, 5),
        cv(eid, 1) + 0.5_ft * cv(eid, 3), cv(eid, 4));
}

// dG2+ (8 in cell, 3 on edge)
template <> inline Eigen::Matrix<FloatType, 1, 3> leftedgeofcell(const DGVector<8>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5_ft * cv(eid, 1) + 1. / 6. * cv(eid, 4),
        cv(eid, 2) - 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 6), cv(eid, 5) - 0.5_ft * cv(eid, 7));
}
template <> inline Eigen::Matrix<FloatType, 1, 3> rightedgeofcell(const DGVector<8>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5_ft * cv(eid, 1) + 1. / 6. * cv(eid, 4),
        cv(eid, 2) + 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 6), cv(eid, 5) + 0.5_ft * cv(eid, 7));
}
template <>
inline Eigen::Matrix<FloatType, 1, 3> bottomedgeofcell(const DGVector<8>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5_ft * cv(eid, 2) + 1. / 6. * cv(eid, 5),
        cv(eid, 1) - 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 7), cv(eid, 4) - 0.5_ft * cv(eid, 6));
}
template <> inline Eigen::Matrix<FloatType, 1, 3> topedgeofcell(const DGVector<8>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5_ft * cv(eid, 2) + 1. / 6. * cv(eid, 5),
        cv(eid, 1) + 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 7), cv(eid, 4) + 0.5_ft * cv(eid, 6));
}

// dQ2 (9 in cell, 3 on edge)
template <> inline Eigen::Matrix<FloatType, 1, 3> leftedgeofcell(const DGVector<9>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5_ft * cv(eid, 1) + 1. / 6. * cv(eid, 4),
        cv(eid, 2) - 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 6),
        cv(eid, 5) - 0.5_ft * cv(eid, 7) + 1. / 6. * cv(eid, 8));
}
template <> inline Eigen::Matrix<FloatType, 1, 3> rightedgeofcell(const DGVector<9>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5_ft * cv(eid, 1) + 1. / 6. * cv(eid, 4),
        cv(eid, 2) + 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 6),
        cv(eid, 5) + 0.5_ft * cv(eid, 7) + 1. / 6. * cv(eid, 8));
}
template <>
inline Eigen::Matrix<FloatType, 1, 3> bottomedgeofcell(const DGVector<9>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) - 0.5_ft * cv(eid, 2) + 1. / 6. * cv(eid, 5),
        cv(eid, 1) - 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 7),
        cv(eid, 4) - 0.5_ft * cv(eid, 6) + 1. / 6. * cv(eid, 8));
}
template <> inline Eigen::Matrix<FloatType, 1, 3> topedgeofcell(const DGVector<9>& cv, size_t eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(cv(eid, 0) + 0.5_ft * cv(eid, 2) + 1. / 6. * cv(eid, 5),
        cv(eid, 1) + 0.5_ft * cv(eid, 3) + 1. / 6. * cv(eid, 7),
        cv(eid, 4) + 0.5_ft * cv(eid, 6) + 1. / 6. * cv(eid, 8));
}

} /* namespace Nextsim */

#endif /* __DGEDGETOOLS_HPP */
