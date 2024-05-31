/*!
 * @file Tools.hpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __TOOLS_HPP
#define __TOOLS_HPP

#include "ParametricTools.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

/*!
 * This namespace collects the auxiliary routines
 */
namespace Tools {

    inline constexpr double SQR(double x) { return x * x; }

    //! Integrates over a DG Vector
    template <int DG> double MeanValue(const ParametricMesh& smesh, const DGVector<DG>& phi)
    {
#define NGP 3
        double mass = 0;
        if (smesh.CoordinateSystem == SPHERICAL) {
            for (int i = 0; i < smesh.nelements; ++i) {
                const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> cos_lat
                    = (ParametricTools::getGaussPointsInElement<NGP>(smesh, i).row(1).array())
                          .cos();
                mass += ((phi.row(i) * PSI<DG, NGP>).array() * cos_lat.array()
                    * (ParametricTools::J<NGP>(smesh, i).array() * GAUSSWEIGHTS<NGP>.array()))
                            .sum();
            }
        } else
            for (int i = 0; i < smesh.nelements; ++i)
                mass += ((phi.row(i) * PSI<DG, NGP>).array()
                    * (ParametricTools::J<NGP>(smesh, i).array() * GAUSSWEIGHTS<NGP>.array()))
                            .sum();

#undef NGP
        return mass;
    }

    static DGVector<1> Landmask(const ParametricMesh& smesh)
    {
        DGVector<1> lm(smesh);
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i)
            lm(i) = smesh.landmask[i] ? 1 : 0;

        return lm;
    }

#define S2A(Q) (Q == 1 ? 1 : (Q == 3 ? 3 : (Q == 6 ? 6 : (Q == 8 ? 6 : -1))))
    template <int DGs>
    DGVector<S2A(DGs)> Delta(const ParametricMesh& smesh, const DGVector<DGs>& E11,
        const DGVector<DGs>& E12, const DGVector<DGs>& E22, const double DeltaMin)
    {
        DGVector<S2A(DGs)> DELTA(smesh);

#define NGP 3

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {
            if (smesh.landmask[i] == 0)
                DELTA.row(i).setZero();
            else {

                const LocalEdgeVector<NGP* NGP> e11_gauss = E11.row(i) * PSI<DGs, NGP>;
                const LocalEdgeVector<NGP* NGP> e12_gauss = E12.row(i) * PSI<DGs, NGP>;
                const LocalEdgeVector<NGP* NGP> e22_gauss = E22.row(i) * PSI<DGs, NGP>;

                const LocalEdgeVector<NGP* NGP> delta_gauss
                    = (DeltaMin * DeltaMin
                          + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                          + 1.50 * e11_gauss.array() * e22_gauss.array()
                          + e12_gauss.array().square())
                          .sqrt()
                          .log10()
                    * ParametricTools::J<NGP>(smesh, i).array() * GAUSSWEIGHTS<NGP>.array();
                DELTA.row(i) = ParametricTools::massMatrix<S2A(DGs)>(smesh, i).inverse()
                    * (PSI<S2A(DGs), NGP> * delta_gauss.transpose());
            }
        }

#undef NGP

        return DELTA;
    }
    template <int DGs>
    DGVector<S2A(DGs)> Shear(const ParametricMesh& smesh, const DGVector<DGs>& E11,
        const DGVector<DGs>& E12, const DGVector<DGs>& E22)
    {
        DGVector<S2A(DGs)> SHEAR(smesh);

#define NGP 3
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {
            if (smesh.landmask[i] == 0)
                SHEAR.row(i).setZero();
            else {

                const LocalEdgeVector<NGP* NGP> e11_gauss = E11.row(i) * PSI<DGs, NGP>;
                const LocalEdgeVector<NGP* NGP> e12_gauss = E12.row(i) * PSI<DGs, NGP>;
                const LocalEdgeVector<NGP* NGP> e22_gauss = E22.row(i) * PSI<DGs, NGP>;

                SHEAR.row(i) = ParametricTools::massMatrix<S2A(DGs)>(smesh, i).inverse()
                    * (PSI<S2A(DGs), NGP>
                        * (((e11_gauss.array() - e22_gauss.array()).square()
                               + 4.0 * e12_gauss.array().square() + 1.e-20)
                                .sqrt()
                                .log10()
                            * ParametricTools::J<NGP>(smesh, i).array() * GAUSSWEIGHTS<NGP>.array())
                              .matrix()
                              .transpose());
            }
        }

#undef NGP

        return SHEAR;
    }

#undef S2A

    template <int DGstress>
    DGVector<DGstress> TensorInvI(const ParametricMesh& smesh, const DGVector<DGstress>& E11,
        const DGVector<DGstress>& E12, const DGVector<DGstress>& E22)
    {
        DGVector<DGstress> Invariant(smesh);

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {
            Invariant.row(i) = 0.5 * (E11.row(i) + E22.row(i));
        }

        return Invariant;
    }

    template <int DGstress>
    DGVector<DGstress> TensorInvII(const ParametricMesh& smesh, const DGVector<DGstress>& E11,
        const DGVector<DGstress>& E12, const DGVector<DGstress>& E22)
    {
        DGVector<DGstress> Invariant(smesh);

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {
            Invariant.row(i)
                = ((0.5 * (E11.row(i) - E22.row(i))).array().pow(2) + E12.row(i).array().pow(2))
                      .sqrt();
        }
        return Invariant;
    }

} /* namespace Tools */

} /* namespace Nextsim */

#endif /* __TOOLS_HPP */
