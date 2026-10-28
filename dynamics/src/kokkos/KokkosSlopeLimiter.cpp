/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosSlopeLimiter.hpp"

namespace Nextsim {

template <int DG>
KokkosSlopeLimiter<DG>::KokkosSlopeLimiter(
    const ParametricMesh& _mesh, const KokkosMesh& kokkosMesh)
    : mesh(kokkosMesh)
{
    CGVector<1> temp;
    temp.resize_by_mesh(_mesh);
    DGVector<1> tempDG;
    tempDG.resize_by_mesh(_mesh);
    if constexpr (DG >= 3) {
        minV = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("minV", temp);
        maxV = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("maxV", temp);
        alpha = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("alpha", tempDG);
        alphaX = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("alphaX", tempDG);
    }
    if constexpr (DG >= 6) {
        dxminV = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("dxminV", temp);
        dxmaxV = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("dxmaxV", temp);
        dyminV = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("dyminV", temp);
        dymaxV = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("dymaxV", temp);
        alphaY = makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>("alphaY", tempDG);
    }
    if constexpr (DG == 8) {
        std::cerr << "No limiting for DG8" << std::endl;
        abort();
    }
}

/*************************************************************/
template <int DG>
void KokkosSlopeLimiter<DG>::initMinMax(DeviceViewCG1& minV, DeviceViewCG1& maxV,
    const ConstDeviceViewDG& phi, DeviceIndex nx, DeviceIndex ny, DeviceIndex comp)
{
    // relative indices of the four vertices in minV/maxV
    const Kokkos::Array<DeviceIndex, 4> cgIndices = { 0, 1, nx + 1, nx + 2 };

    // CPU version uses +-1.e9 but it should probably be min/max
    Kokkos::deep_copy(minV, std::numeric_limits<FloatType>::max());
    Kokkos::deep_copy(maxV, std::numeric_limits<FloatType>::min());

    Kokkos::parallel_for(
        "minMax", phi.extent(0), KOKKOS_LAMBDA(const DeviceIndex eid) {
            const DeviceIndex cx = eid % nx;
            const DeviceIndex cy = eid / nx;
            const DeviceIndex cgi = (nx + 1) * cy + cx;

            const FloatType meanVal = phi(eid, comp);
            for (DeviceIndex j = 0; j < 4; ++j) {
                const DeviceIndex destIdx = cgi + cgIndices[j];
                Kokkos::atomic_min(&minV(destIdx), meanVal);
                Kokkos::atomic_max(&maxV(destIdx), meanVal);
            }
        });
}

/*************************************************************/
// truncates the averages by min or max value
template <int DG> void KokkosSlopeLimiter<DG>::limitMax(const DeviceViewDG& phi, FloatType max)
{
    Kokkos::parallel_for(
        "limitMax", phi.extent(0),
        KOKKOS_LAMBDA(const DeviceIndex c) { phi(c, 0) = Kokkos::min(max, phi(c, 0)); });
}

// truncates the averages by min or max value
template <int DG> void KokkosSlopeLimiter<DG>::limitMin(const DeviceViewDG& phi, FloatType min)
{
    Kokkos::parallel_for(
        "limitMin", phi.extent(0),
        KOKKOS_LAMBDA(const DeviceIndex c) { phi(c, 0) = Kokkos::max(min, phi(c, 0)); });
}

/*************************************************************/
using ConstDeviceViewCG1 = ConstKokkosDeviceView<CGVector<1>>;

KOKKOS_IMPL_FUNCTION static FloatType computeLimit(FloatType midValue,
    const Kokkos::Array<FloatType, 4>& vertexValues, const ConstDeviceViewCG1& _min,
    const ConstDeviceViewCG1& _max, const Kokkos::Array<DeviceIndex, 4>& cgIndices, DeviceIndex cgi)
{
    FloatType al = 1.0; // the limiter
    for (DeviceIndex i = 0; i < 4; ++i) {
        const FloatType dv = vertexValues[i] - midValue; // distance to midpoint
        if (dv > 1.e-8) {
            assert(_max(cgi + cgIndices[i]) >= midValue);
            al = Kokkos::min(al, Kokkos::min(1.0_ft, (_max(cgi + cgIndices[i]) - midValue) / dv));
        }
        if (dv < -1.e-8) {
            assert(_min(cgi + cgIndices[i]) <= midValue);
            al = Kokkos::min(al, Kokkos::min(1.0_ft, (_min(cgi + cgIndices[i]) - midValue) / dv));
        }
        assert(al >= 0);
    }

    return al;
}

template <int DG>
void KokkosSlopeLimiter<DG>::computeAlphas(const DeviceViewDG1& alpha, const ConstDeviceViewDG& phi,
    const ConstDeviceViewCG1& _min, const ConstDeviceViewCG1& _max, DeviceIndex nx, DeviceIndex ny)
{
    // relative indices of the four vertices in minV/maxV
    const Kokkos::Array<DeviceIndex, 4> cgIndices = { 0, 1, nx + 1, nx + 2 };

    assert(alpha.extent(0) == nx * ny);
    Kokkos::parallel_for(
        "computeAlphas", alpha.extent(0), KOKKOS_LAMBDA(const DeviceIndex c) {
            const DeviceIndex cx = c % nx;
            const DeviceIndex cy = c / nx;
            const DeviceIndex cgi = cy * (nx + 1) + cx; // index of lower-left vertex

            // prevent type promotions
            constexpr FloatType off = 0.5;
            // values of phi in the 4 nodes: lower-left, lower-right, upper-left, upper-right
            const Kokkos::Array<FloatType, 4> vertexValues
                = { phi(c, 0) - off * phi(c, 1) - off * phi(c, 2),
                      phi(c, 0) + off * phi(c, 1) - off * phi(c, 2),
                      phi(c, 0) - off * phi(c, 1) + off * phi(c, 2),
                      phi(c, 0) + off * phi(c, 1) + off * phi(c, 2) };

            // value of phi in the midpoint
            const FloatType midValue = phi(c, 0);
            alpha(c) = computeLimit(midValue, vertexValues, _min, _max, cgIndices, cgi);
        });
}

template <int DG>
void KokkosSlopeLimiter<DG>::computeAlphasX(const DeviceViewDG1& alphaX,
    const ConstDeviceViewDG& phi, const ConstDeviceViewCG1& _min, const ConstDeviceViewCG1& _max,
    DeviceIndex nx, DeviceIndex ny)
{
    // relative indices of the four vertices in minV/maxV
    const Kokkos::Array<DeviceIndex, 4> cgIndices = { 0, 1, nx + 1, nx + 2 };

    assert(alphaX.extent(0) == nx * ny);
    Kokkos::parallel_for(
        "computeAlphasX", alphaX.extent(0), KOKKOS_LAMBDA(const DeviceIndex c) {
            const DeviceIndex cx = c % nx;
            const DeviceIndex cy = c / nx;
            const DeviceIndex cgi = cy * (nx + 1) + cx; // index of lower-left vertex

            // prevent type promotions
            constexpr FloatType off = 0.5;
            // values of (d/dx phi) in the 4 nodes: lower-left, lower-right, upper-left, upper-right
            const Kokkos::Array<FloatType, 4> vertexValues = {
                phi(c, 1) - phi(c, 3) - off * phi(c, 5), phi(c, 1) + phi(c, 3) - off * phi(c, 5),
                phi(c, 1) - phi(c, 3) + off * phi(c, 5), phi(c, 1) + phi(c, 3) + off * phi(c, 5)
            };
            // value of d/dx phi in the midpoint
            const FloatType midValue = phi(c, 1);
            alphaX(c) = computeLimit(midValue, vertexValues, _min, _max, cgIndices, cgi);
        });
}

template <int DG>
void KokkosSlopeLimiter<DG>::computeAlphasY(const DeviceViewDG1& alphaY,
    const ConstDeviceViewDG& phi, const ConstDeviceViewCG1& _min, const ConstDeviceViewCG1& _max,
    DeviceIndex nx, DeviceIndex ny)
{
    // relative indices of the four vertices in minV/maxV
    const Kokkos::Array<DeviceIndex, 4> cgIndices = { 0, 1, nx + 1, nx + 2 };

    assert(alphaY.extent(0) == nx * ny);
    Kokkos::parallel_for(
        "computeAlphasY", alphaY.extent(0), KOKKOS_LAMBDA(const DeviceIndex c) {
            const DeviceIndex cx = c % nx;
            const DeviceIndex cy = c / nx;
            const DeviceIndex cgi = cy * (nx + 1) + cx; // index of lower-left vertex

            // prevent type promotions
            constexpr FloatType off = 0.5;
            // values of (d/dx phi) in the 4 nodes: lower-left, lower-right, upper-left, upper-right
            const Kokkos::Array<FloatType, 4> vertexValues = {
                phi(c, 2) - phi(c, 4) - off * phi(c, 5), phi(c, 2) - phi(c, 4) + off * phi(c, 5),
                phi(c, 2) + phi(c, 4) - off * phi(c, 5), phi(c, 2) + phi(c, 4) + off * phi(c, 5)
            };
            // value of d/dx phi in the midpoint
            const FloatType midValue = phi(c, 2);
            alphaY(c) = computeLimit(midValue, vertexValues, _min, _max, cgIndices, cgi);
        });
}

template <int DG>
void KokkosSlopeLimiter<DG>::limitAlphas(
    const DeviceViewDG1& alpha, const DeviceViewDG1& alphaX, const DeviceViewDG1& alphaY)
{
    Kokkos::parallel_for(
        "limitAlphas", alpha.extent(0), KOKKOS_LAMBDA(const DeviceIndex c) {
            alphaX(c) = Kokkos::min(alphaX(c), alphaY(c));
            alpha(c) = Kokkos::max(alpha(c), alphaX(c));
        });
}

template <int DG>
void KokkosSlopeLimiter<DG>::limitHigherOrder(
    const DeviceViewDG& phi, const DeviceViewDG1& alpha, const DeviceViewDG1& alphaX)
{
    Kokkos::parallel_for(
        "limitHigherOrder", phi.extent(0), KOKKOS_LAMBDA(const DeviceIndex c) {
            const FloatType a = alpha(c);
            for (int d = 1; d < 3; ++d)
                phi(c, d) *= a;
            const FloatType aX = alphaX(c);
            for (int d = 3; d < DG; ++d)
                phi(c, d) *= aX;
        });
}

template <int DG> void KokkosSlopeLimiter<DG>::limit(const DeviceViewDG& phi)
{
    if constexpr (DG == 1) // no limiting for dG0
        return;

    // zero order terms & first derivative
    initMinMax(minV, maxV, phi, mesh.nx, mesh.ny); // get max/min values in vertices
    computeAlphas(alpha, phi, minV, maxV, mesh.nx, mesh.ny);

    // derivative & second
    if constexpr (DG == 6) {
        initMinMax(dxminV, dxmaxV, phi, mesh.nx, mesh.ny, 1); // get max/min values in vertices
        initMinMax(dyminV, dymaxV, phi, mesh.nx, mesh.ny, 2); // get max/min values in vertices
        computeAlphasX(alphaX, phi, dxminV, dxmaxV, mesh.nx, mesh.ny);
        computeAlphasY(alphaY, phi, dyminV, dymaxV, mesh.nx, mesh.ny);
        limitAlphas(alpha, alphaX, alphaY);
    }

    limitHigherOrder(phi, alpha, alphaX);
}

template class KokkosSlopeLimiter<DGCOMP>;

}
