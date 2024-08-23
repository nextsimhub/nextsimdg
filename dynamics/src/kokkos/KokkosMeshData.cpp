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

    landMaskDevice = makeKokkosDeviceBitset(mesh.landmask);
}

} // namespace nextsim
