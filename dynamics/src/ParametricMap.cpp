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

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(DG)> dxT
            = ParametricTools::dxT<GAUSSPOINTS1D(DG)>(smesh, eid).array().rowwise()
            * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array();
        const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(DG)> dyT
            = ParametricTools::dyT<GAUSSPOINTS1D(DG)>(smesh, eid).array().rowwise()
            * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array();

        // [J dT^{-T} nabla phi]_1
        AdvectionCellTermX[eid] = PSIx<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * dyT.row(1).array()
            - PSIy<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * dxT.row(1).array();
        // [J dT^{-T} nabla phi]_2

        //! the lat-direction must be scaled with the metric term if in the spherical system
        if (smesh.CoordinateSystem == SPHERICAL) {
            const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(DG)> cos_lat
                = (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(smesh, eid)
                        .row(1)
                        .array())
                      .cos();
            AdvectionCellTermY[eid]
                = PSIy<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * (dxT.row(0).array())
                - PSIx<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * (dyT.row(0).array());
        } else if (smesh.CoordinateSystem == CARTESIAN)
            AdvectionCellTermY[eid]
                = PSIy<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * dxT.row(0).array()
                - PSIx<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * dyT.row(0).array();
        else
            abort();
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
            InverseDGMassMatrix[eid]
                = SphericalTools::massMatrix<DG>(smesh, eid).inverse() / Nextsim::EarthRadius;
    } else if (smesh.CoordinateSystem == CARTESIAN) {
#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid)
            InverseDGMassMatrix[eid] = ParametricTools::massMatrix<DG>(smesh, eid).inverse();
    } else {
        std::cerr << "Coordinate System " << smesh.CoordinateSystem << " not known!" << std::endl;
        abort();
    }
}

//////////////////////////////////////////////////
// Momentum
//////////////////////////////////////////////////

//!
template <int CG> void ParametricMomentumMap<CG>::InitializeLumpedCGMassMatrix()
{
    lumpedcgmass.resize_by_mesh(smesh);

    lumpedcgmass.zero();

#define CGGP(CG) ((CG == 1 ? 1 : 4))

    for (size_t p = 0; p < 2; ++p) // for parallelization
    {
#pragma omp parallel for
        for (size_t iy = p; iy < smesh.ny; iy += 2)
            for (size_t ix = 0; ix < smesh.nx; ++ix) {
                size_t eid = smesh.nx * iy + ix;

                Eigen::Vector<Nextsim::FloatType, (CG == 1 ? 4 : 9)> Meid;

                if (smesh.CoordinateSystem == CARTESIAN) {
                    const Eigen::Matrix<Nextsim::FloatType, 1, CGGP(CG) * CGGP(CG)> J
                        = ParametricTools::J<CGGP(CG)>(smesh, eid).array()
                        * GAUSSWEIGHTS<CGGP(CG)>.array();

                    Meid = PHI<CG, CGGP(CG)> * J.transpose();
                    const Eigen::Matrix<Nextsim::FloatType, 4, 2> coordinates
                        = smesh.coordinatesOfElement(eid);
                } else if (smesh.CoordinateSystem == SPHERICAL) {
                    const Eigen::Matrix<Nextsim::FloatType, 1, CGGP(CG) * CGGP(CG)> cos_lat
                        = (ParametricTools::getGaussPointsInElement<CGGP(CG)>(smesh, eid)
                                .row(1)
                                .array())
                              .cos();

                    const Eigen::Matrix<Nextsim::FloatType, 1, CGGP(CG) * CGGP(CG)> J
                        = ParametricTools::J<CGGP(CG)>(smesh, eid).array()
                        * GAUSSWEIGHTS<CGGP(CG)>.array() * cos_lat.array();

                    Meid = PHI<CG, CGGP(CG)> * J.transpose();
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
}

//!
template <int CG> void ParametricMomentumMap<CG>::InitializeDivSMatrices()
{
    divS1.resize(smesh.nelements);
    divS2.resize(smesh.nelements);
    iMgradX.resize(smesh.nelements);
    iMgradY.resize(smesh.nelements);
    iMJwPSI.resize(smesh.nelements);
    iMJwPSIAdvect.resize(smesh.nelements);
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

        const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(CG2DGSTRESS(CG))> Fx
            = (ParametricTools::dxT<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(smesh, eid).array().rowwise()
                * GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array())
                  .matrix();
        const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(CG2DGSTRESS(CG))> Fy
            = (ParametricTools::dyT<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(smesh, eid).array().rowwise()
                * GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array())
                  .matrix();

        //               [  Fy2  -Fx2 ]
        // A = JF^{-T} = [            ]
        //               [ -Fy1   Fx1 ]
        //

        // the transformed gradient of the CG basis function in the gauss points (OBSERVE SIGN IN
        // SECOND, EIGEN CAN'T START WITH A MINUS
        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), GAUSSPOINTS(CG2DGSTRESS(CG))>
            dx_cg2 = PHIx<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * Fy.row(1).array()
            - PHIy<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * Fx.row(1).array();

        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), GAUSSPOINTS(CG2DGSTRESS(CG))>
            dy_cg2 = PHIy<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * Fx.row(0).array()
            - PHIx<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * Fy.row(0).array();

        // PSI83 is the DG-Basis-function in the guass-point
        // PSI83_{iq} = PSI_i(q)   [ should be called DG_in_GAUSS ]

        const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(CG2DGSTRESS(CG))> J
            = ParametricTools::J<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(smesh, eid);

        iMJwPSIAdvect[eid] = ParametricTools::massMatrix<DGAdvect>(smesh, eid).inverse()
            * (PSI<DGAdvect, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise()
                * (GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array() * J.array()))
                  .matrix();

        if (smesh.CoordinateSystem == CARTESIAN) {
            // divS is used for update of stress (S, nabla Phi) in Momentum
            divS1[eid] = dx_cg2 * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose();
            divS2[eid] = dy_cg2 * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose();

            // iMgradX/Y (inverse-Mass-gradient X/Y) is used to project strain rate from CG to DG
            const Eigen::Matrix<Nextsim::FloatType, CG2DGSTRESS(CG), CG2DGSTRESS(CG)> imass
                = ParametricTools::massMatrix<CG2DGSTRESS(CG)>(smesh, eid).inverse();
            iMgradX[eid] = imass * divS1[eid].transpose();
            iMgradY[eid] = imass * divS2[eid].transpose();

            // imJwPSI is used to compute nonlinear stress update???
            iMJwPSI[eid] = imass
                * (PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise()
                    * (GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array() * J.array()))
                      .matrix();
        } else if (smesh.CoordinateSystem == SPHERICAL) {
            // In spherical coordinates (x,y) coordinates are (lon,lat) coordinates

            const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(CG2DGSTRESS(CG))> cos_lat
                = (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(
                       smesh, eid)
                        .row(1)
                        .array())
                      .cos();
            const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(CG2DGSTRESS(CG))> sin_lat
                = (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(
                       smesh, eid)
                        .row(1)
                        .array())
                      .sin();

            // 1 is lon-derivative, 2 is lat-derivative of the test function
            divS1[eid] = dx_cg2 * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose()
                / Nextsim::EarthRadius;
            divS2[eid] = (dy_cg2.array().rowwise() * cos_lat.array()).matrix()
                * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose()
                / Nextsim::EarthRadius;

            divM[eid] = (PHI<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise()
                            * (J.array() * sin_lat.array()
                                * GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array()))
                            .matrix()
                * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose()
                / Nextsim::EarthRadius;

            const Eigen::Matrix<Nextsim::FloatType, CG2DGSTRESS(CG), CG2DGSTRESS(CG)> imass
                = SphericalTools::massMatrix<CG2DGSTRESS(CG)>(smesh, eid).inverse();
            iMgradX[eid] = imass * divS1[eid].transpose();
            iMgradY[eid] = imass * divS2[eid].transpose();
            iMM[eid] = imass * divM[eid].transpose();

            iMJwPSI[eid] = imass
                * (PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise()
                    * (GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array() * J.array()))
                      .matrix();

        } else
            abort();
    }
}

template class ParametricTransportMap<1>;
template class ParametricTransportMap<3>;
template class ParametricTransportMap<6>;
template class ParametricTransportMap<8>;

template class ParametricMomentumMap<1>;
template class ParametricMomentumMap<2>;

}
