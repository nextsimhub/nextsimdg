/*!
 * @file parametricTools.cpp
 * @date July 28, 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "ParametricTools.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"
#include <Eigen/StdVector>

namespace Nextsim {

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// HOW MANY QUADRATURE POINTS???  I use 3^2 = 9
template <int CG, int DG>
void ParametricTransformation<CG, DG>::BasicInit(const ParametricMesh& smesh)
{
    // resize vectors
    divS1.resize(smesh.nelements);
    divS2.resize(smesh.nelements);
    iMgradX.resize(smesh.nelements);
    iMgradY.resize(smesh.nelements);
    iMJwPSI.resize(smesh.nelements);

    // parallel loop over all elements for computing entries
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {

        //               [  dyT2  -dxT2 ]
        // A = JF^{-T} = [              ]
        //               [ -dyT1   dxT1 ]
        //
        const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(DG)> dxT = (ParametricTools::dxT<GAUSSPOINTS1D(DG)>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array()).matrix();
        const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(DG)> dyT = (ParametricTools::dyT<GAUSSPOINTS1D(DG)>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array()).matrix();

        // the transformed gradient of the CG basis function in the gauss points
        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), GAUSSPOINTS(DG)> dx_cg2 = PHIx<CG, GAUSSPOINTS1D(DG)>.array().rowwise() * dyT.row(1).array() - PHIy<CG, GAUSSPOINTS1D(DG)>.array().rowwise() * dxT.row(1).array();

        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), GAUSSPOINTS(DG)> dy_cg2 = PHIy<CG, GAUSSPOINTS1D(DG)>.array().rowwise() * dxT.row(0).array() - PHIx<CG, GAUSSPOINTS1D(DG)>.array().rowwise() * dyT.row(0).array();

        // PSI83 is the DG-Basis-function in the guass-point
        // PSI83_{iq} = PSI_i(q)   [ should be called DG_in_GAUSS ]

        divS1[eid] = dx_cg2 * PSI<DG, GAUSSPOINTS1D(DG)>.transpose();
        divS2[eid] = dy_cg2 * PSI<DG, GAUSSPOINTS1D(DG)>.transpose();

        const Eigen::Matrix<Nextsim::FloatType, DG, DG> imass = ParametricTools::massMatrix<DG>(smesh, eid).inverse();

        iMgradX[eid] = imass * divS1[eid].transpose();
        iMgradY[eid] = imass * divS2[eid].transpose();

        const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(DG)> J = ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, eid);
        iMJwPSI[eid] = imass * (PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * (GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array() * J.array())).matrix();
    }
}

namespace ParametricTools {
    /*!
     * computes and fills the Q1/Q2 lumped mass matrix
     */
    template <>
    void lumpedCGMassMatrix(const ParametricMesh& smesh,
        CGVector<1>& lumpedcgmass)
    {
        lumpedcgmass.resize_by_mesh(smesh);

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nnodes; ++i)
            lumpedcgmass(i, 0) = 0;

        // parallelization not critical. called just once.

        const size_t sy = smesh.nx + 1;
        size_t i = 0;
        for (size_t iy = 0; iy < smesh.ny; ++iy)
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++i) {
                const double a = smesh.area(i);

                const int n0 = sy * iy + ix;
                lumpedcgmass(n0, 0) += 0.25 * a;
                lumpedcgmass(n0 + 1, 0) += 0.25 * a;
                lumpedcgmass(n0 + sy, 0) += 0.25 * a;
                lumpedcgmass(n0 + sy + 1, 0) += 0.25 * a;
            }
    }

    template <>
    void lumpedCGMassMatrix(const ParametricMesh& smesh,
        CGVector<2>& lumpedcgmass)
    {
        lumpedcgmass.resize_by_mesh(smesh);

        for (size_t i = 0; i < smesh.nnodes; ++i)
            lumpedcgmass(i, 0) = 0;

        // parallelization not critical. called just once.

        const int sy = 2.0 * smesh.nx + 1;
        int i = 0;
        for (size_t iy = 0; iy < smesh.ny; ++iy)
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++i) {
                const double a = smesh.area(i);
                const size_t n0 = 2 * sy * iy + 2 * ix;
                lumpedcgmass(n0, 0) += a / 36.;
                lumpedcgmass(n0 + 2, 0) += a / 36.;
                lumpedcgmass(n0 + 2 * sy, 0) += a / 36.;
                lumpedcgmass(n0 + 2 * sy + 2, 0) += a / 36.;

                lumpedcgmass(n0 + 1, 0) += a / 9.;
                lumpedcgmass(n0 + sy, 0) += a / 9.;
                lumpedcgmass(n0 + sy + 2, 0) += a / 9.;
                lumpedcgmass(n0 + 2 * sy + 1, 0) += a / 9.;

                lumpedcgmass(n0 + sy + 1, 0) += a * 4. / 9.;
            }
    }
}

template class ParametricTransformation<2, 3>;
template class ParametricTransformation<2, 8>;
template class ParametricTransformation<1, 3>;
template class ParametricTransformation<1, 8>;
}
