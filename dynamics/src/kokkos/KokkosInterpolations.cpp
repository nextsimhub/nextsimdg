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
        // much easier to do the pre-computation on CPU
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

    /*************************************************************/
    template <int CG, int DG>
    KokkosDG2CGInterpolator<CG, DG>::KokkosDG2CGInterpolator(const ParametricMesh& smesh)
        : nx(smesh.nx)
        , ny(smesh.ny)
        , PSILagrangeDGCG(PSILagrange<DG, CG + 1>)
    {
    }

    template <int CG, int DG>
    void KokkosDG2CGInterpolator<CG, DG>::operator()(
        const KokkosDeviceView<CGVector<CG>>& dest, const ConstKokkosDeviceView<DGVector<DG>>& src) const
    {
        assert(src.extent(0) == static_cast<int>(nx * ny));
        assert(dest.extent(0) == static_cast<int>((CG * nx + 1) * (CG * ny + 1)));

        Kokkos::deep_copy(dest, 0);

        const DeviceIndex cGDofsPerRow = CG * nx + 1;

        // since all data is needed by the kernel we just capture this
        Kokkos::parallel_for(
            "DG2CG", nx * ny, KOKKOS_CLASS_LAMBDA(const DeviceIndex dgi) {
                const DeviceIndex cy = dgi / nx; //!< y-index of element
                const DeviceIndex cx = dgi % nx; //!< x-index of element
                const DeviceIndex c = dgi;

                const DeviceIndex cgi = CG * cGDofsPerRow * cy + CG * cx;

                auto A = makeEigenMap(src);
                constexpr DeviceIndex NP = (CG + 1) * (CG + 1);
                Eigen::Matrix<FloatType, 1, NP> At = A.row(c) * PSILagrangeDGCG;

                // boundaries
                // top
                if (cy == 0) {
                    for (DeviceIndex i = 0; i < CG + 1; ++i)
                        At(i) *= 2.0;
                }
                // bot
                if (cy == ny - 1) {
                    for (DeviceIndex i = NP - CG - 1; i < NP; ++i)
                        At(i) *= 2.0;
                }
                // left
                if (cx == 0) {
                    for (DeviceIndex i = 0; i < NP; i += CG + 1)
                        At(i) *= 2.0;
                }
                // right
                if (cx == nx - 1) {
                    for (DeviceIndex i = CG; i < NP; i += CG + 1)
                        At(i) *= 2.0;
                }

                // implicit capture outside of
                auto& cg = dest;
                if constexpr (CG == 1) {
                    Kokkos::atomic_add(&cg(cgi), 0.25 * At(0));
                    Kokkos::atomic_add(&cg(cgi + 1), 0.25 * At(1));
                    Kokkos::atomic_add(&cg(cgi + cGDofsPerRow), 0.25 * At(2));
                    Kokkos::atomic_add(&cg(cgi + cGDofsPerRow + 1), 0.25 * At(3));
                } else {
                    Kokkos::atomic_add(&cg(cgi), 0.25 * At(0));
                    Kokkos::atomic_add(&cg(cgi + 1), 0.5 * At(1));
                    Kokkos::atomic_add(&cg(cgi + 2), 0.25 * At(2));
                    Kokkos::atomic_add(&cg(cgi + cGDofsPerRow), 0.5 * At(3));
                    Kokkos::atomic_add(&cg(cgi + cGDofsPerRow + 1), At(4));
                    Kokkos::atomic_add(&cg(cgi + cGDofsPerRow + 2), 0.5 * At(5));
                    Kokkos::atomic_add(&cg(cgi + 2 * cGDofsPerRow), 0.25 * At(6));
                    Kokkos::atomic_add(&cg(cgi + 2 * cGDofsPerRow + 1), 0.5 * At(7));
                    Kokkos::atomic_add(&cg(cgi + 2 * cGDofsPerRow + 2), 0.25 * At(8));
                }
            });
    }

    template class KokkosDG2CGInterpolator<1, 1>;
    template class KokkosDG2CGInterpolator<1, 3>;
    template class KokkosDG2CGInterpolator<1, 6>;
    template class KokkosDG2CGInterpolator<1, 8>;
    template class KokkosDG2CGInterpolator<2, 1>;
    template class KokkosDG2CGInterpolator<2, 3>;
    template class KokkosDG2CGInterpolator<2, 6>;
    template class KokkosDG2CGInterpolator<2, 8>;
}
}