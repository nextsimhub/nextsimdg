/*!
 * @file KokkosParametricMesh.hpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSPARAMETRICMESH_HPP
#define KOKKOSPARAMETRICMESH_HPP

#include "../../include/ParametricMesh.hpp"
#include "KokkosUtils.hpp"

namespace Nextsim {

struct KokkosMeshData {
    KokkosMeshData(const ParametricMesh& mesh);

    // Eigen::Matrix<FloatType, Eigen::Dynamic, 2>
    using ConstDeviceViewVertex = ConstKokkosDeviceView<decltype(ParametricMesh::vertices)>;
    ConstDeviceViewVertex vertices;
    
    using DirichletData = std::array<
        KokkosDeviceMapView<typename decltype(ParametricMesh::dirichlet)::value_type::value_type>,
        ParametricMesh::N_EDGE>;
    DirichletData dirichletDevice;

    ConstDeviceBitset landMaskDevice;

    COORDINATES coordinateSystem;

    KOKKOS_IMPL_DEVICE_FUNCTION Eigen::Matrix<FloatType, 1, 2> edgeVector(DeviceIndex n1, DeviceIndex n2);
};

} // namespace nextsim

#endif // KOKKOSPARAMETRICMESH_HPP