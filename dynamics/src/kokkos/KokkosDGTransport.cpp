/*!
 * @file KokkosDGTransport.cpp
 * @date September 12, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosDGTransport.hpp"

namespace Nextsim{

//! returns the localization of the cell vector to the edges
/*!
 * writes the cell-basis on the edges in the edge basis
 *
 * CELL:
 * DGdegree 0:      1
 * DGdegree 1-2:    + (x-1/2), (y-1/2),
 * DGdegree 3-5:    + (x-1/2)^2-1/12, (y-1/2)^2-1/12, (x-1/2)(y-1/2)
 * DGdegree 6-7:    + (y-1/2)(x-1/2)^2-1/12, (x-1/2)(y-1/2)^2-1/12
 *
 * EDGE:
 * DGdegree 0:    1
 * DGdegree 1: +  (t-1/2)
 * DGdegree 2: +  (t-1/2)^2-1/12
 *
 */
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> leftEdgeOfCell(
    const DGVector<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> rightEdgeOfCell(
    const DGVector<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> bottomEdgeOfCell(
    const DGVector<DG>& cv, DeviceIndex eid);
template <int DG>
Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> topEdgeOfCell(
    const DGVector<DG>& cv, DeviceIndex eid);

namespace Details{
template<typename DGVec>
struct DGDegree;

template<typename T, int DG, class... Properties>
struct DGDegree<Kokkos::View<T*[DG], Properties...>>{
	static constexpr int value = DG;
};
}

template <typename DGVec>
auto leftEdgeOfCell(const DGVec& cv, DeviceIndex eid)
{
	constexpr int DG = Details::DGDegree<DGVec>::value;
	if constexpr (DG == 1) {
    	return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
	} else {
		return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
	}
}

// dG0 (1 in cell, 1 on edge)
/*
Eigen::Matrix<FloatType, 1, 1> leftEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}*/

Eigen::Matrix<FloatType, 1, 1> rightEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}

Eigen::Matrix<FloatType, 1, 1> bottomEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}
 Eigen::Matrix<FloatType, 1, 1> topEdgeOfCell(const DGVector<1>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 1>(cv(eid, 0));
}

// dG1 (3 in cell, 2 on edge)

Eigen::Matrix<FloatType, 1, 2> leftEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}

Eigen::Matrix<FloatType, 1, 2> rightEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}

Eigen::Matrix<FloatType, 1, 2> bottomEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
 Eigen::Matrix<FloatType, 1, 2> topEdgeOfCell(const DGVector<3>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)

Eigen::Matrix<FloatType, 1, 3> leftEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3), cv(eid, 2) - 0.5 * cv(eid, 5),
        cv(eid, 4));
}

Eigen::Matrix<FloatType, 1, 3> rightEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3), cv(eid, 2) + 0.5 * cv(eid, 5),
        cv(eid, 4));
}

Eigen::Matrix<FloatType, 1, 3> bottomEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4), cv(eid, 1) - 0.5 * cv(eid, 5),
        cv(eid, 3));
}
 Eigen::Matrix<FloatType, 1, 3> topEdgeOfCell(const DGVector<6>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4), cv(eid, 1) + 0.5 * cv(eid, 5),
        cv(eid, 3));
}

// dG2+ (8 in cell, 3 on edge)

Eigen::Matrix<FloatType, 1, 3> leftEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) - 0.5 * cv(eid, 7));
}

Eigen::Matrix<FloatType, 1, 3> rightEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6), cv(eid, 4) + 0.5 * cv(eid, 7));
}

Eigen::Matrix<FloatType, 1, 3> bottomEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) - 0.5 * cv(eid, 6));
}
 Eigen::Matrix<FloatType, 1, 3> topEdgeOfCell(const DGVector<8>& cv, DeviceIndex eid)
{
    return Eigen::Matrix<FloatType, 1, 3>(
        cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), cv(eid, 3) + 0.5 * cv(eid, 6));
}

}