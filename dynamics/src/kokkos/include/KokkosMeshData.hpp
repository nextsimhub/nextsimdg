/*!
 * @file KokkosParametricMesh.hpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSPARAMETRICMESH_HPP
#define KOKKOSPARAMETRICMESH_HPP

#include "../../include/ParametricMesh.hpp"
#include "KokkosUtils.hpp"

#include <Kokkos_Bitset.hpp>

namespace Nextsim {

struct KokkosMeshData {
    KokkosMeshData(const ParametricMesh& mesh);

    std::array<KokkosDeviceMapView<size_t>, 4> dirichletDevice;
    Kokkos::ConstBitset<Kokkos::DefaultExecutionSpace> landMaskDevice;
};

} // namespace nextsim

#endif // KOKKOSPARAMETRICMESH_HPP