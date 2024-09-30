/*!
 * @file    KokkosDGTransport.hpp
 * @date    September 26 2024
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSINTERPOLATIONS_HPP
#define __KOKKOSINTERPOLATIONS_HPP

#include "../include/cgVector.hpp"
#include "../include/codeGenerationDGinGauss.hpp"
#include "../include/dgVector.hpp"
#include "KokkosUtils.hpp"
// #include "../include/ParametricMesh.hpp"

namespace Nextsim {
namespace Interpolations {

    // This functor keeps its own precomputed map that could be very large.
    // Therefore you should use only one instance of this.
    template <int DG, int CG> class KokkosCG2DGInterpolator {
    public:
        KokkosCG2DGInterpolator(const ParametricMesh& smesh);

        void operator()(KokkosDeviceView<DGVector<DG>>& dgDevice,
            const ConstKokkosDeviceView<CGVector<CG>>& cgDevice) const;

    private:
        using CG2DGMatrix = Eigen::Matrix<FloatType, DG, CG == 2 ? 9 : 4>;
        KokkosDeviceMapView<CG2DGMatrix> cG2DGMatrixDevice;
        DeviceIndex nx;
        DeviceIndex ny;
        DeviceIndex nelements;
    };

    template <int CG, int DG> class KokkosDG2CGInterpolator {
    public:
        KokkosDG2CGInterpolator(const ParametricMesh& smesh);

        void operator()(const KokkosDeviceView<CGVector<CG>>& dest,
            const ConstKokkosDeviceView<DGVector<DG>>& src) const;

    private:
        using PSILagrangeDGCGType = decltype(PSILagrange<DG, CG + 1>);
        PSILagrangeDGCGType PSILagrangeDGCG;
        DeviceIndex nx;
        DeviceIndex ny;
    };
}
}

#endif // __KOKKOSINTERPOLATIONS_HPP