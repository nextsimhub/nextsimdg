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
template <>
void ParametricTransformation<2, 8>::BasicInit(const SasipMesh& smesh)
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
        const Eigen::Matrix<Nextsim::FloatType, 2, 9> dxT = (ParametricTools::dxT<3>(smesh, eid).array().rowwise() * GAUSSWEIGHTS_3.array()).matrix();
        const Eigen::Matrix<Nextsim::FloatType, 2, 9> dyT = (ParametricTools::dyT<3>(smesh, eid).array().rowwise() * GAUSSWEIGHTS_3.array()).matrix();

        // the transformed gradient of the CG basis function in the gauss points
        const Eigen::Matrix<Nextsim::FloatType, 9, 9> dx_cg2 = CG_CG2_dx_in_GAUSS3.array().rowwise() * dyT.row(1).array() - CG_CG2_dy_in_GAUSS3.array().rowwise() * dxT.row(1).array();

        const Eigen::Matrix<Nextsim::FloatType, 9, 9> dy_cg2 = CG_CG2_dy_in_GAUSS3.array().rowwise() * dxT.row(0).array() - CG_CG2_dx_in_GAUSS3.array().rowwise() * dyT.row(0).array();

        // PSI83 is the DG-Basis-function in the guass-point
        // PSI83_{iq} = PSI_i(q)   [ should be called DG_in_GAUSS ]

        divS1[eid] = dx_cg2 * PSI<8,3>.transpose();
        divS2[eid] = dy_cg2 * PSI<8,3>.transpose();

        const Eigen::Matrix<Nextsim::FloatType, 8, 8> imass = ParametricTools::massMatrix<8>(smesh, eid).inverse();

        iMgradX[eid] = imass * divS1[eid].transpose();
        iMgradY[eid] = imass * divS2[eid].transpose();

        const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, eid);
        iMJwPSI[eid] = imass * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS_3.array() * J.array())).matrix();
    }
}

namespace ParametricTools {
    /*!
     * computes and fills the Q1/Q2 lumped mass matrix
     */
    template <>
    void lumpedCGMassMatrix(const SasipMesh& smesh,
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
    void lumpedCGMassMatrix(const SasipMesh& smesh,
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

template class ParametricTransformation<2, 8>;
// template class ParametricTransformation<1, 3>;
}
