/*!
 * @file KokkosMeshData.cpp
 * @date August 22, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosMeshData.hpp"

namespace Nextsim {

KokkosMeshData::KokkosMeshData(const ParametricMesh& mesh)
{
    // boundary data
    for (int i = 0; i < ParametricMesh::Edge::N_EDGE; ++i) {
        dirichletDevice[i]
            = makeKokkosDeviceViewMap("dirichlet" + std::to_string(i), mesh.dirichlet[i], true);
    }

    // unfortunately there is no more direct way to initialize a Kokkos::Bitset
    const unsigned nBits = mesh.landmask.size();
    Kokkos::Bitset<Kokkos::HostSpace> landMaskHost(nBits);
    landMaskHost.clear();
    for (unsigned i = 0; i < nBits; ++i) {
        if (mesh.landmask[i])
            landMaskHost.set(i);
    }
    // mutable device bitset -> const device bitset
    Kokkos::Bitset<Kokkos::DefaultExecutionSpace> landMaskTemp(nBits);
    Kokkos::deep_copy(landMaskTemp, landMaskHost);
    landMaskDevice = landMaskTemp;
}

} // namespace nextsim
