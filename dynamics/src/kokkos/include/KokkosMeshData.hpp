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

    using DirichletData = std::array<
        KokkosDeviceMapView<typename decltype(ParametricMesh::dirichlet)::value_type::value_type>,
        ParametricMesh::N_EDGE>;
    DirichletData dirichletDevice;
    ConstDeviceBitset landMaskDevice;
};

} // namespace nextsim

#endif // KOKKOSPARAMETRICMESH_HPP