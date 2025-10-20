/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSSLOPELIMITER_HPP
#define __KOKKOSSLOPELIMITER_HPP

#include "KokkosMesh.hpp"
#include "KokkosUtils.hpp"
#include "ParametricMesh.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"

namespace Nextsim {

// Slope limiter based on D. Kuzmin: A vertex-based hierarchical slope limiter for p-adaptive
// discontinuous Galerkin methods. Journal of Computational and Applied Mathematics 2010.

template <int DG> class KokkosSlopeLimiter {
    const KokkosMesh& mesh;

    using DeviceViewCG1 = KokkosDeviceView<CGVector<1>>;
    using ConstDeviceViewCG1 = ConstKokkosDeviceView<CGVector<1>>;
    using DeviceViewDG1 = KokkosDeviceView<DGVector<1>>;
    using DeviceViewDG = KokkosDeviceView<DGVector<DG>>;
    using ConstDeviceViewDG = ConstKokkosDeviceView<DGVector<DG>>;

    //! minimum and maximum values at the mesh nodes
    DeviceViewCG1 minV, maxV;
    DeviceViewCG1 dxminV, dxmaxV;
    DeviceViewCG1 dyminV, dymaxV;

    //! alpha-values for limiting
    DeviceViewDG1 alpha, alphaX, alphaY;

public:
    KokkosSlopeLimiter(const ParametricMesh& _mesh, const KokkosMesh& kokkosMesh);

    // performs the vertex-based limiting
    void limit(const DeviceViewDG& phi);

    static void limitMax(const DeviceViewDG& phi, FloatType max);
    static void limitMin(const DeviceViewDG& phi, FloatType min);

    // ************************************************************* //
    // implementation (treat as protected)
    // gets minimum and maximum values at all mesh vertices
    static void initMinMax(DeviceViewCG1& minV, DeviceViewCG1& maxV, const ConstDeviceViewDG& phi,
        DeviceIndex nx, DeviceIndex ny, DeviceIndex comp = 0);

    static void computeAlphas(const DeviceViewDG1& alpha, const ConstDeviceViewDG& phi,
        const ConstDeviceViewCG1& _min, const ConstDeviceViewCG1& _max, DeviceIndex nx,
        DeviceIndex ny);

    static void computeAlphasX(const DeviceViewDG1& alphaX, const ConstDeviceViewDG& phi,
        const ConstDeviceViewCG1& _min, const ConstDeviceViewCG1& _max, DeviceIndex nx,
        DeviceIndex ny);

    static void computeAlphasY(const DeviceViewDG1& alphaY, const ConstDeviceViewDG& phi,
        const ConstDeviceViewCG1& _min, const ConstDeviceViewCG1& _max, DeviceIndex nx,
        DeviceIndex ny);

    static void limitAlphas(
        const DeviceViewDG1& alpha, const DeviceViewDG1& alphaX, const DeviceViewDG1& alphaY);
    static void limitHigherOrder(
        const DeviceViewDG& phi, const DeviceViewDG1& alpha, const DeviceViewDG1& alphaX);
};

}

#endif
