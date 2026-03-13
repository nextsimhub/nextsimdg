/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosMesh.hpp"

namespace Nextsim {

KokkosMesh::KokkosMesh(const ParametricMesh& mesh)
    : coordinateSystem(mesh.CoordinateSystem)
    , nx(mesh.nx)
    , ny(mesh.ny)
{
    // vertices
    vertices = makeKokkosDeviceView("vertices", mesh.vertices, MakeViewOptions::DEVICE_COPY);

    landMaskDevice = makeKokkosDeviceBitset(mesh.landmask);
}

} // namespace nextsim
