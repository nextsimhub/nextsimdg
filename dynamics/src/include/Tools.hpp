/*!
 * @file Tools.hpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __TOOLS_HPP
#define __TOOLS_HPP

#include "dgVector.hpp"

namespace Nextsim {

/*!
 * This namespace collects the auxiliary routines
 */
namespace Tools {

    inline constexpr double SQR(double x) { return x * x; }

    template <int DGstress>
    CellVector<1> Delta(const Mesh& mesh, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const double DeltaMin)
    {
        CellVector<1> DELTA(mesh);

#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            DELTA(i, 0) = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
        }

        return DELTA;
    }

    template <int DGstress>
    CellVector<1> Shear(const Mesh& mesh, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
			const CellVector<DGstress>& E22)
    {
      CellVector<1> SHEAR(mesh);

#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            SHEAR(i, 0) = sqrt((SQR(E11(i, 0) - E22(i, 0)) + 4.0 * SQR(E12(i, 0))));
        }

        return SHEAR;
    }

    template <int DGstress>
    CellVector<1> Delta(const SasipMesh& smesh, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const double DeltaMin)
    {
        CellVector<1> DELTA(smesh);

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            DELTA(i, 0) = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
        }

        return DELTA;
    }

    template <int DGstress>
    CellVector<1> Shear(const SasipMesh& smesh, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
			const CellVector<DGstress>& E22)
    {
        CellVector<1> SHEAR(smesh);

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {
            SHEAR(i, 0) = sqrt((SQR(E11(i, 0) - E22(i, 0)) + 4.0 * SQR(E12(i, 0))));
        }

        return SHEAR;
    }

    template <int DGstress, int DGtracer>
    void TensorInvariants(const Mesh& mesh, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, CellVector<DGtracer>& Inv1, CellVector<DGtracer>& Inv2)
    {

#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            // \dot E_I
            Inv1(i, 0) = 0.5 * (E11(i, 0) + E22(i, 0));
            // \dot E_II
            Inv2(i, 0) = std::hypot(0.5 * (E11(i, 0) - E22(i, 0)), E12(i, 0));
        }
    }

    template <int DGstress>
    void TensorInvariants(const Mesh& mesh, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, CellVector<DGstress>& Inv1, CellVector<DGstress>& Inv2)
    {

#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            // \dot E_I
            Inv1.row(i) = 0.5 * (E11.row(i) + E22.row(i));
            // \dot E_II
            Inv2.row(i) = ((0.5 * (E11.row(i) - E22.row(i))).array().pow(2) + E12.row(i).array().pow(2)).sqrt();
        }
    }

    template <int DGstress, int DGtracer>
    void ElastoParams(const Mesh& mesh, const CellVector<DGstress>& E11,
        const CellVector<DGstress>& E12, const CellVector<DGstress>& E22,
        const CellVector<DGtracer>& H, const CellVector<DGtracer>& A, const double DeltaMin,
        const double Pstar, CellVector<DGtracer>& MU1, CellVector<DGtracer>& MU2)
    {

#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            double DELTA = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            if (DELTA <= 0)
                std::cout << DELTA << std::endl;

            assert(DELTA >= 0);

            //! Ice strength
            double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

            double zeta = P / 2.0 / DELTA;
            double eta = zeta / 4;

            MU1(i, 0) = 2 * eta;
            MU2(i, 0) = zeta - eta;
        }
    }

} /* namespace Tools */

} /* namespace Nextsim */

#endif /* __TOOLS_HPP */
