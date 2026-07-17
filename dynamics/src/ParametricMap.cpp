/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#include "ParametricMap.hpp"
#include "ParametricTools.hpp"
#include "VectorManipulations.hpp"

namespace Nextsim {

template <int DG> void ParametricTransportMap<DG>::InitializeAdvectionCellTerms()
{
    // for advection
    AdvectionCellTermX.clear();
    AdvectionCellTermY.clear();

    AdvectionCellTermX.resize(smesh.nelements);
    AdvectionCellTermY.resize(smesh.nelements);

    // gradient of transformation
    //      [ dxT1, dyT1 ]     //            [ dyT2, -dxT2 ]
    // dT = 		     // J dT^{-T}=
    //      [ dxT2, dyT2 ]     //            [ -dyT1, dxT1 ]
    //
    // given as [dxT1, dxT2, dyT1, dyT2] ->  [dyT2, -dxT2, -dyT1, dxT1 ]

    // J dT^{-T} nabla Phi  = [dyT2 * PSIx - dxT2 * PSIy, -dyT1 * PSIx + dxT1 * PSIy]
    // PSIx, PSIy are DG x QQ - matrices
    // dxT, dyT are 2 x QQ - matrices

    // Store wq * phi(q)

    constexpr int GP_DG = GAUSSPOINTS(DG);
    constexpr int GP_1D = GAUSSPOINTS1D(DG);

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        const Eigen::Matrix<double, 2, GP_DG> dxT
            = ParametricTools::dxT<GP_1D, double>(smesh, eid).array().rowwise()
            * GAUSSWEIGHTS<GP_1D, double>.array();
        const Eigen::Matrix<double, 2, GP_DG> dyT
            = ParametricTools::dyT<GP_1D, double>(smesh, eid).array().rowwise()
            * GAUSSWEIGHTS<GP_1D, double>.array();

        // [J dT^{-T} nabla phi]_1
        AdvectionCellTermX[eid] = (PSIx<DG, GP_1D, double>.array().rowwise() * dyT.row(1).array()
            - PSIy<DG, GP_1D, double>.array().rowwise() * dxT.row(1).array())
                                      .template cast<FloatType>();
        // [J dT^{-T} nabla phi]_2

        //! the lat-direction must be scaled with the metric term if in the spherical system
        if (smesh.CoordinateSystem == SPHERICAL) {
            const Eigen::Matrix<double, 1, GP_DG> cos_lat
                = (ParametricTools::getGaussPointsInElement<GP_1D, double>(smesh, eid)
                        .row(1)
                        .array())
                      .cos();
            AdvectionCellTermY[eid]
                = (PSIy<DG, GP_1D, double>.array().rowwise() * (dxT.row(0).array())
                    - PSIx<DG, GP_1D, double>.array().rowwise() * (dyT.row(0).array()))
                      .template cast<FloatType>();
        } else if (smesh.CoordinateSystem == CARTESIAN) {
            AdvectionCellTermY[eid]
                = (PSIy<DG, GP_1D, double>.array().rowwise() * dxT.row(0).array()
                    - PSIx<DG, GP_1D, double>.array().rowwise() * dyT.row(0).array())
                      .template cast<FloatType>();
        } else {
            abort();
        }
    }
}

template <int DG> void ParametricTransportMap<DG>::InitializeInverseDGMassMatrix()
{
    // for advection
    InverseDGMassMatrix.clear();
    InverseDGMassMatrix.resize(smesh.nelements);

    if (smesh.CoordinateSystem == SPHERICAL) {
#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid)
            InverseDGMassMatrix[eid] = (SphericalTools::massMatrix<DG, double>(smesh, eid).inverse()
                / Nextsim::EarthRadius)
                                           .template cast<FloatType>();
    } else if (smesh.CoordinateSystem == CARTESIAN) {
#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid)
            InverseDGMassMatrix[eid] = ParametricTools::massMatrix<DG, double>(smesh, eid)
                                           .inverse()
                                           .template cast<FloatType>();
    } else {
        std::cerr << "Coordinate System " << smesh.CoordinateSystem << " not known!" << std::endl;
        abort();
    }
}

//////////////////////////////////////////////////
// Momentum
//////////////////////////////////////////////////

//!
template <int CG, int DG> void ParametricMomentumMap<CG, DG>::InitializeLumpedCGMassMatrix()
{
    // Compute lumped mass matric for cG(CG)

    lumpedcgmass.resize_by_mesh(smesh);

    lumpedcgmass.zero();

#define CGGP(CG) ((CG == 1 ? 1 : 4))

    for (size_t p = 0; p < 2; ++p) // for parallelization
    {
#pragma omp parallel for
        for (size_t iy = p; iy < smesh.ny; iy += 2)
            for (size_t ix = 0; ix < smesh.nx; ++ix) {
                size_t eid = smesh.nx * iy + ix;

                Eigen::Vector<double, (CG == 1 ? 4 : 9)> Meid;

                if (smesh.CoordinateSystem == CARTESIAN) {
                    const Eigen::Matrix<double, 1, CGGP(CG) * CGGP(CG)> J
                        = ParametricTools::J<CGGP(CG), double>(smesh, eid).array()
                        * GAUSSWEIGHTS<CGGP(CG), double>.array();

                    Meid = PHI<CG, CGGP(CG), double> * J.transpose();
                } else if (smesh.CoordinateSystem == SPHERICAL) {
                    const Eigen::Matrix<double, 1, CGGP(CG) * CGGP(CG)> cos_lat
                        = (ParametricTools::getGaussPointsInElement<CGGP(CG), double>(smesh, eid)
                                .row(1)
                                .array())
                              .cos();

                    const Eigen::Matrix<double, 1, CGGP(CG) * CGGP(CG)> J
                        = ParametricTools::J<CGGP(CG), double>(smesh, eid).array()
                        * GAUSSWEIGHTS<CGGP(CG), double>.array() * cos_lat.array();

                    Meid = PHI<CG, CGGP(CG), double> * J.transpose();
                } else
                    abort();

                // index of first dof in element
                const size_t sy = CG * smesh.nx + 1;
                const size_t n0 = CG * iy * sy + CG * ix;

                if (CG == 1) {
                    lumpedcgmass(n0, 0) += Meid(0);
                    lumpedcgmass(n0 + 1, 0) += Meid(1);

                    lumpedcgmass(n0 + sy, 0) += Meid(2);
                    lumpedcgmass(n0 + 1 + sy, 0) += Meid(3);
                } else if (CG == 2) {
                    lumpedcgmass(n0, 0) += Meid(0);
                    lumpedcgmass(n0 + 1, 0) += Meid(1);
                    lumpedcgmass(n0 + 2, 0) += Meid(2);

                    lumpedcgmass(n0 + sy, 0) += Meid(3);
                    lumpedcgmass(n0 + 1 + sy, 0) += Meid(4);
                    lumpedcgmass(n0 + 2 + sy, 0) += Meid(5);

                    lumpedcgmass(n0 + 2 * sy, 0) += Meid(6);
                    lumpedcgmass(n0 + 1 + 2 * sy, 0) += Meid(7);
                    lumpedcgmass(n0 + 2 + 2 * sy, 0) += Meid(8);
                } else
                    abort();
            }
    }
    //    VectorManipulations::CGAddPeriodic(smesh, lumpedcgmass);

    // Build the cG1 mass matrix. If CG=1 this is redundant. But CG=2 is standard.
    lumpedcg1mass.resize_by_mesh(smesh);
    lumpedcg1mass.zero();

    for (size_t p = 0; p < 2; ++p) // for parallelization
    {
#pragma omp parallel for
        for (size_t iy = p; iy < smesh.ny; iy += 2)
            for (size_t ix = 0; ix < smesh.nx; ++ix) {
                size_t eid = smesh.nx * iy + ix;

                Eigen::Vector<double, 4> Meid;

                if (smesh.CoordinateSystem == CARTESIAN) {
                    const Eigen::Matrix<double, 1, 2 * 2> J
                        = ParametricTools::J<2, double>(smesh, eid).array()
                        * GAUSSWEIGHTS<2, double>.array();

                    Meid = PHI<1, 2, double> * J.transpose();
                } else if (smesh.CoordinateSystem == SPHERICAL) {
                    const Eigen::Matrix<double, 1, 2 * 2> cos_lat
                        = (ParametricTools::getGaussPointsInElement<2, double>(smesh, eid)
                                .row(1)
                                .array())
                              .cos();

                    const Eigen::Matrix<double, 1, 2 * 2> J
                        = ParametricTools::J<2, double>(smesh, eid).array()
                        * GAUSSWEIGHTS<2, double>.array() * cos_lat.array();

                    Meid = PHI<1, 2, double> * J.transpose();
                } else
                    abort();

                // index of first dof in element
                const size_t sy = smesh.nx + 1;
                const size_t n0 = iy * sy + ix;

                lumpedcg1mass(n0, 0) += Meid(0);
                lumpedcg1mass(n0 + 1, 0) += Meid(1);

                lumpedcg1mass(n0 + sy, 0) += Meid(2);
                lumpedcg1mass(n0 + 1 + sy, 0) += Meid(3);
            }
    }
}

//!
template <int CG, int DG> void ParametricMomentumMap<CG, DG>::InitializeDivSMatrices()
{
    dX_SSH.resize(smesh.nelements);
    dY_SSH.resize(smesh.nelements);
    divS1.resize(smesh.nelements);
    divS2.resize(smesh.nelements);
    iMgradX.resize(smesh.nelements);
    iMgradY.resize(smesh.nelements);
    iMJwPSI.resize(smesh.nelements);
    iMJwPSI_dam.resize(smesh.nelements);
    if (smesh.CoordinateSystem == SPHERICAL) {
        divM.resize(smesh.nelements);
        iMM.resize(smesh.nelements);
    }

    // parallel loop over all elements for computing entries
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {

        //     [ Fx   Fx ]
        // F = [         ]
        //     [ Fy   Fy ]

        constexpr int GP = GAUSSPOINTS(CG2DGSTRESS(CG));
        constexpr int GP_1D = GAUSSPOINTS1D(CG2DGSTRESS(CG));
        constexpr int DG_STRESS = CG2DGSTRESS(CG);

        const Eigen::Matrix<double, 2, GP> Fx
            = (ParametricTools::dxT<GP_1D, double>(smesh, eid).array().rowwise()
                * GAUSSWEIGHTS<GP_1D, double>.array())
                  .matrix();
        const Eigen::Matrix<double, 2, GP> Fy
            = (ParametricTools::dyT<GP_1D, double>(smesh, eid).array().rowwise()
                * GAUSSWEIGHTS<GP_1D, double>.array())
                  .matrix();

        //               [  Fy2  -Fx2 ]
        // A = JF^{-T} = [            ]
        //               [ -Fy1   Fx1 ]
        //

        // the transformed gradient of the CG basis function in the gauss points (OBSERVE SIGN IN
        // SECOND, EIGEN CAN'T START WITH A MINUS
        const Eigen::Matrix<double, (CG == 2 ? 9 : 4), GP> dx_cg2
            = PHIx<CG, GP_1D, double>.array().rowwise() * Fy.row(1).array()
            - PHIy<CG, GP_1D, double>.array().rowwise() * Fx.row(1).array();

        const Eigen::Matrix<double, (CG == 2 ? 9 : 4), GP> dy_cg2
            = PHIy<CG, GP_1D, double>.array().rowwise() * Fx.row(0).array()
            - PHIx<CG, GP_1D, double>.array().rowwise() * Fy.row(0).array();

        // same but using CG1 basis functions. Required for seaSurfaceHeight
        const Eigen::Matrix<double, 4, GP> dx_cg1
            = PHIx<1, GP_1D, double>.array().rowwise() * Fy.row(1).array()
            - PHIy<1, GP_1D, double>.array().rowwise() * Fx.row(1).array();

        const Eigen::Matrix<double, 4, GP> dy_cg1
            = PHIy<1, GP_1D, double>.array().rowwise() * Fx.row(0).array()
            - PHIx<1, GP_1D, double>.array().rowwise() * Fy.row(0).array();

        const Eigen::Matrix<double, 1, GP> J = ParametricTools::J<GP_1D, double>(smesh, eid);

        if (smesh.CoordinateSystem == CARTESIAN) {
            // divS is used for update of stress (S, nabla Phi) in Momentum
            const Eigen::Matrix<double, CGDOFS(CG), DG_STRESS> divS1_
                = dx_cg2 * PSI<DG_STRESS, GP_1D, double>.transpose();
            const Eigen::Matrix<double, CGDOFS(CG), DG_STRESS> divS2_
                = dy_cg2 * PSI<DG_STRESS, GP_1D, double>.transpose();

            divS1[eid] = divS1_.template cast<FloatType>();
            divS2[eid] = divS2_.template cast<FloatType>();

            // dX_SSH and dY_SSH are used to compute the gradient of the sea surface height
            // they store (d_[x/y] PHI_j, PHI_i)
            dX_SSH[eid] = (dx_cg1 * PHI<1, GP_1D, double>.transpose()).template cast<FloatType>();
            dY_SSH[eid] = (dy_cg1 * PHI<1, GP_1D, double>.transpose()).template cast<FloatType>();

            // iMgradX/Y (inverse-Mass-gradient X/Y) is used to project strain rate from CG to DG
            const Eigen::Matrix<double, DG_STRESS, DG_STRESS> imass
                = ParametricTools::massMatrix<DG_STRESS, double>(smesh, eid).inverse();
            iMgradX[eid] = (imass * divS1_.transpose()).template cast<FloatType>();
            iMgradY[eid] = (imass * divS2_.transpose()).template cast<FloatType>();

            // imJwPSI is used to compute nonlinear stress update
            iMJwPSI[eid] = (imass
                * (PSI<DG_STRESS, GP_1D, double>.array().rowwise()
                    * (GAUSSWEIGHTS<GP_1D, double>.array() * J.array()))
                      .matrix())
                               .template cast<FloatType>();
            // same but for the damage. However, we use the same number of Gausspoints as
            // for the DG-stress variant above for easier use in BBM Stress update
            const Eigen::Matrix<double, DG, DG> imass_dam
                = ParametricTools::massMatrix<DG, double>(smesh, eid).inverse();
            iMJwPSI_dam[eid] = (imass_dam
                * (PSI<DG, GP_1D, double>.array().rowwise()
                    * (GAUSSWEIGHTS<GP_1D, double>.array() * J.array()))
                      .matrix())
                                   .template cast<FloatType>();

        } else if (smesh.CoordinateSystem == SPHERICAL) {
            // In spherical coordinates (x,y) coordinates are (lon,lat) coordinates

            const Eigen::Matrix<double, 1, GP> cos_lat
                = (ParametricTools::getGaussPointsInElement<GP_1D, double>(smesh, eid)
                        .row(1)
                        .array())
                      .cos();
            const Eigen::Matrix<double, 1, GP> sin_lat
                = (ParametricTools::getGaussPointsInElement<GP_1D, double>(smesh, eid)
                        .row(1)
                        .array())
                      .sin();

            // 1 is lon-derivative, 2 is lat-derivative of the test function
            const Eigen::Matrix<double, CGDOFS(CG), DG_STRESS> divS1_
                = dx_cg2 * PSI<DG_STRESS, GP_1D, double>.transpose() / Nextsim::EarthRadius;
            const Eigen::Matrix<double, CGDOFS(CG), DG_STRESS> divS2_
                = (dy_cg2.array().rowwise() * cos_lat.array()).matrix()
                * PSI<DG_STRESS, GP_1D, double>.transpose() / Nextsim::EarthRadius;

            divS1[eid] = divS1_.template cast<FloatType>();
            divS2[eid] = divS2_.template cast<FloatType>();

            const Eigen::Matrix<double, CGDOFS(CG), DG_STRESS> divM_
                = (PHI<CG, GP_1D, double>.array().rowwise()
                      * (J.array() * sin_lat.array() * GAUSSWEIGHTS<GP_1D, double>.array()))
                      .matrix()
                * PSI<DG_STRESS, GP_1D, double>.transpose() / Nextsim::EarthRadius;

            divM[eid] = divM_.template cast<FloatType>();

            // same for CG1 (Sea-Surface Height)
            dX_SSH[eid] = (dx_cg1 * PHI<1, GP_1D, double>.transpose() / Nextsim::EarthRadius)
                              .template cast<FloatType>();
            dY_SSH[eid] = ((dy_cg1.array().rowwise() * cos_lat.array()).matrix()
                * PHI<1, GP_1D, double>.transpose() / Nextsim::EarthRadius)
                              .template cast<FloatType>();

            const Eigen::Matrix<double, DG_STRESS, DG_STRESS> imass
                = SphericalTools::massMatrix<DG_STRESS, double>(smesh, eid).inverse();
            iMgradX[eid] = (imass * divS1_.transpose()).template cast<FloatType>();
            iMgradY[eid] = (imass * divS2_.transpose()).template cast<FloatType>();
            iMM[eid] = (imass * divM_.transpose()).template cast<FloatType>();

            // same as in Cartesian but using spherical mass matrix and adding cosine
            // for proper scale of the integral
            iMJwPSI[eid] = (imass
                * (PSI<DG_STRESS, GP_1D, double>.array().rowwise()
                    * (GAUSSWEIGHTS<GP_1D, double>.array() * J.array() * cos_lat.array()))
                      .matrix())
                               .template cast<FloatType>();

            // smae for DG advection (damage)
            const Eigen::Matrix<double, DG, DG> imass_dam
                = SphericalTools::massMatrix<DG, double>(smesh, eid).inverse();
            iMJwPSI_dam[eid] = (imass_dam
                * (PSI<DG, GP_1D, double>.array().rowwise()
                    * (GAUSSWEIGHTS<GP_1D, double>.array() * J.array() * cos_lat.array()))
                      .matrix())
                                   .template cast<FloatType>();

        } else
            abort();
    }
}

//!
template <int CG, int DG> void ParametricMomentumMap<CG, DG>::InitializeCGLandmask()
{
    cglandmask.resize_by_mesh(smesh);
    for (size_t i = 0; i < cglandmask.size(); ++i)
        cglandmask[i] = 1.0;

    size_t eid = 0;
    const size_t dofsinrow = CG * smesh.nx + 1;
    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        size_t cgid = iy * CG * dofsinrow; // first cg index in row
        for (size_t ix = 0; ix < smesh.nx; ++ix, ++eid, cgid += CG) {
            if (smesh.landmask[eid] == 0) // on land
                for (size_t cy = 0; cy <= CG; ++cy)
                    for (size_t cx = 0; cx <= CG; ++cx)
                        cglandmask[cgid + cy * dofsinrow + cx] = 0;
        }
    }
}

template class ParametricTransportMap<1>;
template class ParametricTransportMap<3>;
template class ParametricTransportMap<6>;
template class ParametricTransportMap<8>;

template class ParametricMomentumMap<1, 1>;
template class ParametricMomentumMap<2, 1>;
template class ParametricMomentumMap<1, 3>;
template class ParametricMomentumMap<2, 3>;
template class ParametricMomentumMap<1, 6>;
template class ParametricMomentumMap<2, 6>;

}
