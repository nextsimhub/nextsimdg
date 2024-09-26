#include "include/KokkosInterpolations.hpp"
#include "../include/ParametricTools.hpp"

namespace Nextsim {
namespace Interpolations {

    template <int DG, int CG>
    KokkosCG2DGInterpolator<DG, CG>::KokkosCG2DGInterpolator(const ParametricMesh& smesh)
        : nx(smesh.nx)
        , ny(smesh.ny)
        , nelements(smesh.nelements)
    {
        std::vector<CG2DGMatrix> cG2DGMatrix;
        cG2DGMatrix.resize(smesh.nelements);
		#pragma omp parallel for
        for (size_t dgi = 0; dgi < smesh.nelements; ++dgi) {
            if (smesh.CoordinateSystem == CARTESIAN) {
                cG2DGMatrix[dgi]
                    = ((ParametricTools::massMatrix<DG>(smesh, dgi).inverse()
                           * (PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise()
                               * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, dgi).array()
                                   * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array())
                                     .array())
                                 .matrix())
                        * PHI<CG, GAUSSPOINTS1D(DG)>.transpose());
            } else {
                cG2DGMatrix[dgi]
                    = ((ParametricTools::massMatrix<DG>(smesh, dgi).inverse()
                           * (PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise()
                               * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, dgi).array()
                                   * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array()
                                   * ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(
                                       smesh, dgi)
                                         .row(1)
                                         .array()
                                         .cos())
                                     .array())
                                 .matrix())
                        * PHI<CG, GAUSSPOINTS1D(DG)>.transpose());
            }
        }
        cG2DGMatrixDevice = makeKokkosDeviceViewMap("cG2DGMatrix", cG2DGMatrix, true);
    }

    template <int DG, int CG>
    void KokkosCG2DGInterpolator<DG, CG>::operator()(KokkosDeviceView<DGVector<DG>>& dgDevice,
        const ConstKokkosDeviceView<CGVector<CG>>& cgDevice) const
    {
        assert((CG * nx + 1) * (CG * ny + 1) == cgDevice.extent(0));
        assert(nelements == dgDevice.extent(0));

        const int cgshift = CG * nx + 1; //!< Index shift for each row

        // since all data is needed by the kernel we just capture this
        Kokkos::parallel_for(
            "CG2DG", nx * ny, KOKKOS_CLASS_LAMBDA(const DeviceIndex dgi) {
                const DeviceIndex iy = dgi / nx; //!< y-index of element
                const DeviceIndex ix = dgi % nx; //!< x-index of element
                const DeviceIndex cgi
                    = CG * cgshift * iy + CG * ix; //!< lower/left Index in cg vector

                Eigen::Matrix<FloatType, (CG == 2 ? 9 : 4), 1>
                    cgLocal; //!< the local unknowns in the element
                // we need to use cgDevice outside of the constexpr branch for implicit capture
                const auto& cg = cgDevice;
                if constexpr (CG == 1) {
                    cgLocal << cg(cgi), cg(cgi + 1), cg(cgi + cgshift), cg(cgi + 1 + cgshift);
                } else {
                    cgLocal << cg(cgi), cg(cgi + 1), cg(cgi + 2), cg(cgi + cgshift),
                        cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift), cg(cgi + 2 * cgshift),
                        cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift);
                }
                // solve:  (Vdg, PHI) = (Vcg, PHI) with mapping to spher. coord.
                auto dg = makeEigenMap(dgDevice);
                dg.row(dgi) = cG2DGMatrixDevice[dgi] * cgLocal;
            });
    }

    template class KokkosCG2DGInterpolator<1, 1>;
    template class KokkosCG2DGInterpolator<3, 1>;
    template class KokkosCG2DGInterpolator<6, 1>;
    template class KokkosCG2DGInterpolator<8, 1>;
    template class KokkosCG2DGInterpolator<1, 2>;
    template class KokkosCG2DGInterpolator<3, 2>;
    template class KokkosCG2DGInterpolator<6, 2>;
    template class KokkosCG2DGInterpolator<8, 2>;

}
}