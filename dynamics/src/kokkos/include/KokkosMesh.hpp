/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSPARAMETRICMESH_HPP
#define KOKKOSPARAMETRICMESH_HPP

#include "../../../core/src/kokkos/include/KokkosUtils.hpp"
#include "../../include/ParametricMesh.hpp"

namespace Nextsim {

struct KokkosMesh {
    KokkosMesh(const ParametricMesh& mesh);

    using ConstDeviceViewVertex = ConstKokkosDeviceView<decltype(ParametricMesh::vertices)>;
    ConstDeviceViewVertex vertices;

    ConstDeviceBitset landMaskDevice;

    COORDINATES coordinateSystem;
    DeviceIndex nx;
    DeviceIndex ny;

    // has to be defined in the header because of linking issues
    // (could work with relocatable device code but at a performance hit)
    KOKKOS_IMPL_FUNCTION inline Eigen::Matrix<FloatType, 1, 2> edgeVector(
        DeviceIndex n1, DeviceIndex n2) const
    {
        Eigen::Matrix<FloatType, 1, 2> dv;
        dv(0) = vertices(n2, 0) - vertices(n1, 0);
        dv(1) = vertices(n2, 1) - vertices(n1, 1);
        //! In spherical coordinates (and greenland) we must check for the Pi -> -Pi jump
        if (coordinateSystem == SPHERICAL) {
            if (dv(0, 0) > 0.5_ft * M_PI)
                dv(0, 0) -= 2.0_ft * M_PI;
            if (dv(0, 0) < -0.5_ft * M_PI)
                dv(0, 0) += 2.0_ft * M_PI;
        }

        return dv;
    }
};

} // namespace nextsim

#endif // KOKKOSPARAMETRICMESH_HPP
