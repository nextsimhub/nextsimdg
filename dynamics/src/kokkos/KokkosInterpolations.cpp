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
                    = ((SphericalTools::massMatrix<DG>(smesh, dgi).inverse()
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
        cG2DGMatrixDevice
            = makeKokkosDeviceViewMap("cG2DGMatrix", cG2DGMatrix, MakeViewOptions::ALWAYS_COPY);
    }

    template <int DG, int CG>
    void KokkosCG2DGInterpolator<DG, CG>::operator()(const KokkosDeviceView<DGVector<DG>>& dgDevice,
        const ConstKokkosDeviceView<CGVector<CG>>& cgDevice) const
    {
        assert((CG * nx + 1) * (CG * ny + 1) == cgDevice.extent(0));
        assert(nelements == dgDevice.extent(0));

        const DeviceIndex cgshift = CG * nx + 1; //!< Index shift for each row

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
    void KokkosDG2CGInterpolator<CG, DG>::operator()(const KokkosDeviceView<CGVector<CG>>& dest,
        const ConstKokkosDeviceView<DGVector<DG>>& src) const
    {
        assert(src.extent(0) == static_cast<int>(nx * ny));
        assert(dest.extent(0) == static_cast<int>((CG * nx + 1) * (CG * ny + 1)));

        // explicit execution space enables asynchronous execution
        auto execSpace = Kokkos::DefaultExecutionSpace();
        Kokkos::deep_copy(execSpace, dest, 0.0);

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

                // for implicit capture dest has to be used outside of the constexpr branch
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

    /*************************************************************/
    void kokkosCG12CG2(const KokkosDeviceView<CGVector<2>>& dest,
        const ConstKokkosDeviceView<CGVector<1>>& src, const DeviceIndex nx, const DeviceIndex ny)
    {
        assert(src.extent(0) == (nx + 1) * (ny + 1));
        assert(dest.extent(0) == (2 * nx + 1) * (2 * ny + 1));

        const DeviceIndex cg1Row = nx + 1;
        Kokkos::parallel_for(
            "CG12CG2", src.extent(0), KOKKOS_LAMBDA(const DeviceIndex icg1) {
                const DeviceIndex iy = icg1 / cg1Row;
                const DeviceIndex ix = icg1 % cg1Row;
                const DeviceIndex icg2 = (2 * nx + 1) * 2 * iy + 2 * ix;

                // outer nodes
                const FloatType dof0 = src(icg1);
                dest(icg2) = dof0;

                // along horizontal lines
                FloatType dof1;
                if (ix < nx) {
                    dof1 = src(icg1 + 1);
                    dest(icg2 + 1) = 0.5 * (dof0 + dof1);
                }

                // along vertical lines
                FloatType dof2;
                if (iy < ny) {
                    dof2 = src(icg1 + cg1Row);
                    dest(icg2 + 2 * nx + 1) = 0.5 * (dof0 + dof2);
                }

                // midpoints
                if (ix < nx && iy < ny) {
                    const FloatType dof3 = src(icg1 + cg1Row + 1);
                    dest(icg2 + 2 * nx + 1 + 1) = 0.25 * (dof0 + dof1 + dof2 + dof3);
                }
            });
    }
}
}