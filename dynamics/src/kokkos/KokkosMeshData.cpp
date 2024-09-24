/*!
 * @file KokkosMeshData.cpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosMeshData.hpp"

namespace Nextsim {

KokkosMeshData::KokkosMeshData(const ParametricMesh& mesh)
    : coordinateSystem(mesh.CoordinateSystem)
    , nx(mesh.nx)
    , ny(mesh.ny)
{
    // vertices
    vertices = makeKokkosDeviceView("vertices", mesh.vertices, true);

    // boundary data
    for (int i = 0; i < ParametricMesh::Edge::N_EDGE; ++i) {
        dirichletDevice[i]
            = makeKokkosDeviceViewMap("dirichlet" + std::to_string(i), mesh.dirichlet[i], true);
    }

    landMaskDevice = makeKokkosDeviceBitset(mesh.landmask);
}
/*
KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, 1, 2> KokkosMeshData::edgeVector(
    DeviceIndex n1, DeviceIndex n2) const
{
    Eigen::Matrix<FloatType, 1, 2> dv;
    dv(0) = vertices(n2, 0) - vertices(n1, 0);
    dv(1) = vertices(n2, 1) - vertices(n2, 0);
    //    Eigen::Matrix<Nextsim::FloatType, 1, 2> dv
    //        = vertices.block<1, 2>(n2, 0) - vertices.block<1, 2>(n1, 0);
    //! In spherical coordinates (and greenland) we must check for the Pi -> -Pi jump
    if (coordinateSystem == SPHERICAL) {
        if (dv(0, 0) > 0.5 * M_PI)
            dv(0, 0) -= 2.0 * M_PI;
        if (dv(0, 0) < -0.5 * M_PI)
            dv(0, 0) += 2.0 * M_PI;
    }

    return dv;
}*/

} // namespace nextsim
