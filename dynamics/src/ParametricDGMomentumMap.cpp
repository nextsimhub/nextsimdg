/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 *
 * Precomputed matrices of the LDG momentum discretization (ParametricDGMomentumMap)
 */

#include "ParametricMap.hpp"
#include "ParametricTools.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgEdgeTools.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>

namespace Nextsim {

//! aborts unless the mesh is CARTESIAN (the LDG matrices are not derived for SPHERICAL yet)
static void requireCartesian(const ParametricMesh& smesh, const std::string& function)
{
    if (smesh.CoordinateSystem != CARTESIAN) {
        std::cerr << "ParametricDGMomentumMap::" << function
                  << ": only CARTESIAN coordinates are supported" << std::endl;
        abort();
    }
}

template <int DG, int DGstress, int DGstrain>
template <int DGweight>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::InitializeWeightedDGMassMatrix(
    const DGVector<DGweight>& weight)
{
    requireCartesian(smesh, "InitializeWeightedDGMassMatrix");

    WeightedDGMassMatrix.resize(smesh.nelements);
    InverseWeightedDGMassMatrix.resize(smesh.nelements);

    // the integrand PSI_i * PSI_j * weight has degree 2*DGDEGREE(DG) + DGDEGREE(DGweight)
    constexpr int NGP = NGAUSSFORDEGREE1D(2 * DGDEGREE(DG) + DGDEGREE(DGweight));

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        const Eigen::Matrix<FloatType, 1, NGP * NGP> J = ParametricTools::J<NGP>(smesh, eid);
        const Eigen::Matrix<FloatType, 1, NGP * NGP> w_gauss = weight.row(eid) * PSI<DGweight, NGP>;

        WeightedDGMassMatrix[eid] = (PSI<DG, NGP>.array().rowwise()
                                        * (GAUSSWEIGHTS<NGP>.array() * J.array() * w_gauss.array()))
                                        .matrix()
            * PSI<DG, NGP>.transpose();
        InverseWeightedDGMassMatrix[eid] = WeightedDGMassMatrix[eid].inverse();
    }
}

template <int DG, int DGstress, int DGstrain>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::InitializeDGDivMatrices()
{
    requireCartesian(smesh, "InitializeDGDivMatrices");

    dgDivS1.resize(smesh.nelements);
    dgDivS2.resize(smesh.nelements);
    iMgradTauX.resize(smesh.nelements);
    iMgradTauY.resize(smesh.nelements);
    InverseStrainMassMatrix.resize(smesh.nelements);

    // the integrands (d PSI_DG) * PSI_DGstress and (d PSI_DGstrain) * PSI_DG have degree
    // DGDEGREE(DG) + DGDEGREE(DGstress), assuming DGDEGREE(DGstrain) <= DGDEGREE(DGstress)
    static_assert(DGDEGREE(DGstrain) <= DGDEGREE(DGstress));
    constexpr int NGP = NGAUSSFORDEGREE1D(DGDEGREE(DG) + DGDEGREE(DGstress));

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        const Eigen::Matrix<FloatType, 2, NGP * NGP> Fx
            = (ParametricTools::dxT<NGP>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<NGP>.array())
                  .matrix();
        const Eigen::Matrix<FloatType, 2, NGP * NGP> Fy
            = (ParametricTools::dyT<NGP>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<NGP>.array())
                  .matrix();

        // the transformed gradient of the velocity basis in the gauss points (OBSERVE SIGN IN
        // SECOND, EIGEN CAN'T START WITH A MINUS)
        const Eigen::Matrix<FloatType, DG, NGP * NGP> dx_v
            = PSIx<DG, NGP>.array().rowwise() * Fy.row(1).array()
            - PSIy<DG, NGP>.array().rowwise() * Fx.row(1).array();
        const Eigen::Matrix<FloatType, DG, NGP * NGP> dy_v
            = PSIy<DG, NGP>.array().rowwise() * Fx.row(0).array()
            - PSIx<DG, NGP>.array().rowwise() * Fy.row(0).array();

        // the transformed gradient of the strain basis in the gauss points
        const Eigen::Matrix<FloatType, DGstrain, NGP * NGP> dx_tau
            = PSIx<DGstrain, NGP>.array().rowwise() * Fy.row(1).array()
            - PSIy<DGstrain, NGP>.array().rowwise() * Fx.row(1).array();
        const Eigen::Matrix<FloatType, DGstrain, NGP * NGP> dy_tau
            = PSIy<DGstrain, NGP>.array().rowwise() * Fx.row(0).array()
            - PSIx<DGstrain, NGP>.array().rowwise() * Fy.row(0).array();

        // (sigma, nabla v)
        dgDivS1[eid] = dx_v * PSI<DGstress, NGP>.transpose();
        dgDivS2[eid] = dy_v * PSI<DGstress, NGP>.transpose();

        // M_eps^-1 (u, div tau)
        InverseStrainMassMatrix[eid] = ParametricTools::massMatrix<DGstrain>(smesh, eid).inverse();
        iMgradTauX[eid] = InverseStrainMassMatrix[eid] * (dx_tau * PSI<DG, NGP>.transpose());
        iMgradTauY[eid] = InverseStrainMassMatrix[eid] * (dy_tau * PSI<DG, NGP>.transpose());
    }
}

template <int DG, int DGstress, int DGstrain>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::InitializeJwPSI()
{
    requireCartesian(smesh, "InitializeJwPSI");

    JwPSI.resize(smesh.nelements);
    constexpr int NGP = GAUSSPOINTS1D(DG);

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        JwPSI[eid] = PSI<DG, NGP>.array().rowwise()
            * (GAUSSWEIGHTS<NGP>.array() * ParametricTools::J<NGP>(smesh, eid).array());
}

template <int DG, int DGstress, int DGstrain>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::InitializeIMJwPSI()
{
    requireCartesian(smesh, "InitializeIMJwPSI");
    // the inverse stress mass matrix is the inverse strain mass matrix
    static_assert(DGstress == DGstrain);
    assert(InverseStrainMassMatrix.size() == smesh.nelements);

    iMJwPSI.resize(smesh.nelements);
    // same Gauss points as used in the stress update
    constexpr int NGP = GAUSSPOINTS1D(DGstress);

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        iMJwPSI[eid] = InverseStrainMassMatrix[eid]
            * (PSI<DGstress, NGP>.array().rowwise()
                * (GAUSSWEIGHTS<NGP>.array() * ParametricTools::J<NGP>(smesh, eid).array()))
                  .matrix();
}

/*!
 * Values of the cell basis functions in the Gauss points of the side `side` (0 bottom,
 * 1 right, 2 top, 3 left): result(g, j) = phi_j(gauss point g). Computed with the same
 * edge traces and PSIe tables as the flux terms, such that the Gauss points are ordered as
 * in PSIe_w.
 */
template <int DG, int NGP> static Eigen::Matrix<FloatType, NGP, DG> cellBasisOnSide(const int side)
{
    Eigen::Matrix<FloatType, NGP, DG> values;
    DGVector<DG> unitCoefficients;
    unitCoefficients.resize(1, DG); // a "mesh" of one cell

    for (int j = 0; j < DG; ++j) {
        unitCoefficients.setZero();
        unitCoefficients(0, j) = 1.0;

        Eigen::Matrix<FloatType, 1, EDGEDOFS(DG)> edgeCoefficients;
        if (side == 0)
            edgeCoefficients = bottomedgeofcell<DG>(unitCoefficients, 0);
        else if (side == 1)
            edgeCoefficients = rightedgeofcell<DG>(unitCoefficients, 0);
        else if (side == 2)
            edgeCoefficients = topedgeofcell<DG>(unitCoefficients, 0);
        else
            edgeCoefficients = leftedgeofcell<DG>(unitCoefficients, 0);

        values.col(j) = (edgeCoefficients * PSIe<EDGEDOFS(DG), NGP>).transpose();
    }
    return values;
}

//! As cellBasisOnSide, multiplied with the Gauss weights (= PSIe_w)
template <int DG, int NGP>
static Eigen::Matrix<FloatType, NGP, DG> weightedCellBasisOnSide(const int side)
{
    if (side == 0)
        return PSIe_w<DG, NGP, 0>;
    else if (side == 1)
        return PSIe_w<DG, NGP, 1>;
    else if (side == 2)
        return PSIe_w<DG, NGP, 2>;
    else
        return PSIe_w<DG, NGP, 3>;
}

template <int DG, int DGstress, int DGstrain>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::InitializeJumpPenalty(
    const FloatType a, const FloatType b)
{
    requireCartesian(smesh, "InitializeJumpPenalty");
    if (!(a < 0.5)) {
        std::cerr << "ParametricDGMomentumMap::InitializeJumpPenalty: the Dirichlet penalty "
                     "b/(h*(0.5-a)) needs a < 0.5 (got a = "
                  << a << ")" << std::endl;
        abort();
    }

    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;

    // 1. Trace mass matrices on the reference edge, the integrand has degree 2*DGDEGREE(DG)
    constexpr int NGP = NGAUSSFORDEGREE1D(2 * DGDEGREE(DG));
    for (int side = 0; side < 4; ++side) {
        const int oppositeSide = (side + 2) % 4;
        TraceMassSameCell[side]
            = weightedCellBasisOnSide<DG, NGP>(side).transpose() * cellBasisOnSide<DG, NGP>(side);
        TraceMassNeighborCell[side] = weightedCellBasisOnSide<DG, NGP>(side).transpose()
            * cellBasisOnSide<DG, NGP>(oppositeSide);
    }

    // length of a (straight) edge
    auto lengthOfEdgeX = [&](const size_t edge) {
        const Eigen::Matrix<FloatType, 2, 2> ends = smesh.coordinatesOfEdgeX(edge);
        return (ends.row(1) - ends.row(0)).norm();
    };
    auto lengthOfEdgeY = [&](const size_t edge) {
        const Eigen::Matrix<FloatType, 2, 2> ends = smesh.coordinatesOfEdgeY(edge);
        return (ends.row(1) - ends.row(0)).norm();
    };
    auto isOcean = [&](const size_t cell) { return smesh.landmask[cell] != 0; };

    // 2. Interior edges, w_e = b * l_e / (0.5*(h(c1)+h(c2))): coupling block
    //    w_e * TraceMassNeighborCell (test c1, trial c2) and the w_e * TraceMassSameCell
    //    contributions to the diagonal blocks of c1 and c2.
    //    No coupling on boundary edges and edges touching land.
    JumpPenaltyDiagonalBlock.assign(smesh.nelements, Matrix<DG, DG>::Zero());

    JumpPenaltyCouplingX.assign(nx * (ny + 1), Matrix<DG, DG>::Zero());
    for (size_t iy = 0; iy + 1 < ny; ++iy)
        for (size_t ix = 0; ix < nx; ++ix) {
            const size_t cellBelow = iy * nx + ix;
            const size_t cellAbove = cellBelow + nx;
            if (!isOcean(cellBelow) || !isOcean(cellAbove))
                continue;
            const size_t edge = cellAbove; // bottom edge of cellAbove
            const FloatType weight
                = b * lengthOfEdgeX(edge) / (0.5_ft * (smesh.h(cellBelow) + smesh.h(cellAbove)));
            JumpPenaltyCouplingX[edge] = weight * TraceMassNeighborCell[2]; // top of cellBelow
            JumpPenaltyDiagonalBlock[cellBelow] += weight * TraceMassSameCell[2]; // top
            JumpPenaltyDiagonalBlock[cellAbove] += weight * TraceMassSameCell[0]; // bottom
        }

    JumpPenaltyCouplingY.assign((nx + 1) * ny, Matrix<DG, DG>::Zero());
    for (size_t iy = 0; iy < ny; ++iy)
        for (size_t ix = 0; ix + 1 < nx; ++ix) {
            const size_t cellLeft = iy * nx + ix;
            const size_t cellRight = cellLeft + 1;
            if (!isOcean(cellLeft) || !isOcean(cellRight))
                continue;
            const size_t edge = iy * (nx + 1) + ix + 1; // right edge of cellLeft
            const FloatType weight
                = b * lengthOfEdgeY(edge) / (0.5_ft * (smesh.h(cellLeft) + smesh.h(cellRight)));
            JumpPenaltyCouplingY[edge] = weight * TraceMassNeighborCell[1]; // right of cellLeft
            JumpPenaltyDiagonalBlock[cellLeft] += weight * TraceMassSameCell[1]; // right
            JumpPenaltyDiagonalBlock[cellRight] += weight * TraceMassSameCell[3]; // left
        }

    // 3. Diagonal blocks: plus the Dirichlet edges
    addDirichletJumpPenalty(a, b);
}

template <int DG, int DGstress, int DGstrain>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::addDirichletJumpPenalty(
    const FloatType a, const FloatType b)
{
    const size_t nx = smesh.nx;

    // w_e = b * l_e / (h(c) * (0.5-a)).
    // smesh.dirichlet[side] lists the cells having a Dirichlet edge on that side.
    // Sequential: a corner cell is in two lists.
    for (int side = 0; side < 4; ++side)
        for (const size_t cell : smesh.dirichlet[side]) {
            if (smesh.landmask[cell] == 0)
                continue;
            const size_t ix = cell % nx;
            const size_t iy = cell / nx;
            Eigen::Matrix<FloatType, 2, 2> ends;
            if (side == 0)
                ends = smesh.coordinatesOfEdgeX(cell);
            else if (side == 1)
                ends = smesh.coordinatesOfEdgeY(iy * (nx + 1) + ix + 1);
            else if (side == 2)
                ends = smesh.coordinatesOfEdgeX(cell + nx);
            else
                ends = smesh.coordinatesOfEdgeY(iy * (nx + 1) + ix);
            const FloatType length = (ends.row(1) - ends.row(0)).norm();
            const FloatType weight = b * length / (smesh.h(cell) * (0.5_ft - a));
            JumpPenaltyDiagonalBlock[cell] += weight * TraceMassSameCell[side];
        }
}

template <int DG, int DGstress, int DGstrain>
template <int DGadvection>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::UpdateViscousJumpPenalty(const FloatType a,
    const FloatType bMin, const FloatType bMax, const FloatType C, const VPParameters& params,
    const DGVector<DGstrain>& E11, const DGVector<DGstrain>& E12, const DGVector<DGstrain>& E22,
    const DGVector<DGadvection>& H, const DGVector<DGadvection>& A)
{
    assert(JumpPenaltyCouplingX.size() == smesh.nx * (smesh.ny + 1));
    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;

    // zeta is not a polynomial: one more Gauss point than for phi_i * phi_j
    constexpr int NGP = NGAUSSFORDEGREE1D(2 * DGDEGREE(DG) + 2);
    using EdgeValues = Eigen::Matrix<FloatType, NGP, 1>;

    // basis functions in the Gauss points of each side: velocity (also weighted), strain,
    // ice height and concentration
    std::array<Eigen::Matrix<FloatType, NGP, DG>, 4> basis, weightedBasis;
    std::array<Eigen::Matrix<FloatType, NGP, DGstrain>, 4> basisStrain;
    std::array<Eigen::Matrix<FloatType, NGP, DGadvection>, 4> basisAdvection;
    for (int side = 0; side < 4; ++side) {
        basis[side] = cellBasisOnSide<DG, NGP>(side);
        weightedBasis[side] = weightedCellBasisOnSide<DG, NGP>(side);
        basisStrain[side] = cellBasisOnSide<DGstrain, NGP>(side);
        basisAdvection[side] = cellBasisOnSide<DGadvection, NGP>(side);
    }

    // trace of zeta = P / (2 Delta) of `cell` on its side `side`, P and Delta as in
    // mEVPStressUpdate (the viscosities use P, also with the replacement pressure)
    const FloatType sqrDeltaMin = params.deltaMin * params.deltaMin;
    auto zetaOnSide = [&](const size_t cell, const int side) -> EdgeValues {
        const EdgeValues e11 = basisStrain[side] * E11.row(cell).transpose();
        const EdgeValues e12 = basisStrain[side] * E12.row(cell).transpose();
        const EdgeValues e22 = basisStrain[side] * E22.row(cell).transpose();
        const EdgeValues h
            = (basisAdvection[side] * H.row(cell).transpose()).array().max(0.0_ft).matrix();
        const EdgeValues conc = (basisAdvection[side] * A.row(cell).transpose())
                                    .array()
                                    .max(0.0_ft)
                                    .min(1.0_ft)
                                    .matrix();
        const EdgeValues P
            = (params.pStar * h.array() * (params.compactionParam * (1.0_ft - conc.array())).exp())
                  .matrix();
        const EdgeValues DELTA
            = (sqrDeltaMin + 1.25_ft * (e11.array().square() + e22.array().square())
                + 1.50_ft * e11.array() * e22.array() + e12.array().square())
                  .sqrt()
                  .matrix();
        return (0.5_ft * P.array() / DELTA.array()).matrix();
    };

    // int_ref w(xi) phi_i^(sideTest) phi_j^(sideTrial) dxi (w contains the edge length)
    auto edgeBlock = [&](const EdgeValues& w, const int sideTest, const int sideTrial) {
        return (weightedBasis[sideTest].transpose() * w.asDiagonal() * basis[sideTrial]).eval();
    };

    // geometric factor l_e / h_perp of an interior edge, h_perp = min(|c1|,|c2|) / l_e the
    // width of the thinner neighbour normal to the edge (sqrt(area) overestimates it on
    // stretched or distorted cells)
    auto edgeFactor = [&](const Eigen::Matrix<FloatType, 2, 2>& ends, const size_t c1,
                          const size_t c2) {
        return (ends.row(1) - ends.row(0)).squaredNorm() / std::min(smesh.area(c1), smesh.area(c2));
    };

    // min(max(C {zeta}, bMin), bMax) in the Gauss points of an interior edge
    auto cappedViscosity = [&](const EdgeValues& zeta1, const EdgeValues& zeta2) -> EdgeValues {
        return (C * 0.5_ft * (zeta1 + zeta2)).array().max(bMin).min(bMax).matrix();
    };

    auto isOcean = [&](const size_t cell) { return smesh.landmask[cell] != 0; };

    JumpPenaltyDiagonalBlock.assign(smesh.nelements, Matrix<DG, DG>::Zero());

    // interior X-edges, {zeta} from the cell below (its top, side 2) and the cell above (its
    // bottom, side 0). Parallel over the columns: two columns never share a cell.
#pragma omp parallel for
    for (size_t ix = 0; ix < nx; ++ix)
        for (size_t iy = 0; iy + 1 < ny; ++iy) {
            const size_t cellBelow = iy * nx + ix;
            const size_t cellAbove = cellBelow + nx;
            const size_t edge = cellAbove; // bottom edge of cellAbove
            if (!isOcean(cellBelow) || !isOcean(cellAbove))
                continue; // coupling stays 0 (InitializeJumpPenalty)
            const EdgeValues w = edgeFactor(smesh.coordinatesOfEdgeX(edge), cellBelow, cellAbove)
                * cappedViscosity(zetaOnSide(cellBelow, 2), zetaOnSide(cellAbove, 0));
            JumpPenaltyCouplingX[edge] = edgeBlock(w, 2, 0);
            JumpPenaltyDiagonalBlock[cellBelow] += edgeBlock(w, 2, 2);
            JumpPenaltyDiagonalBlock[cellAbove] += edgeBlock(w, 0, 0);
        }

        // interior Y-edges, {zeta} from the cell on the left (its right, side 1) and the cell on
        // the right (its left, side 3). Parallel over the rows.
#pragma omp parallel for
    for (size_t iy = 0; iy < ny; ++iy)
        for (size_t ix = 0; ix + 1 < nx; ++ix) {
            const size_t cellLeft = iy * nx + ix;
            const size_t cellRight = cellLeft + 1;
            const size_t edge = iy * (nx + 1) + ix + 1; // right edge of cellLeft
            if (!isOcean(cellLeft) || !isOcean(cellRight))
                continue;
            const EdgeValues w = edgeFactor(smesh.coordinatesOfEdgeY(edge), cellLeft, cellRight)
                * cappedViscosity(zetaOnSide(cellLeft, 1), zetaOnSide(cellRight, 3));
            JumpPenaltyCouplingY[edge] = edgeBlock(w, 1, 3);
            JumpPenaltyDiagonalBlock[cellLeft] += edgeBlock(w, 1, 1);
            JumpPenaltyDiagonalBlock[cellRight] += edgeBlock(w, 3, 3);
        }

    // Dirichlet edges: constant bMax, not zeta which vanishes in the shear zones along the
    // walls and would no longer enforce u = 0
    addDirichletJumpPenalty(a, bMax);
}

template <int DG, int DGstress, int DGstrain>
void ParametricDGMomentumMap<DG, DGstress, DGstrain>::ApplyMassPlusJumpPenalty(
    const FloatType massCoefficient, const FloatType penaltyCoefficient, const DGVector<DG>& x,
    DGVector<DG>& y) const
{
    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;

#pragma omp parallel for
    for (size_t cell = 0; cell < smesh.nelements; ++cell) {
        const size_t ix = cell % nx;
        const size_t iy = cell / nx;
        const size_t leftEdge = iy * (nx + 1) + ix;

        // J x restricted to this cell: own block ...
        Eigen::Matrix<FloatType, DG, 1> penalty
            = JumpPenaltyDiagonalBlock[cell] * x.row(cell).transpose();

        // ... minus the coupling with each neighbour (zero block = no coupling): the edge
        // block as is when the cell is c1 (below / left of the edge), transposed when it is c2
        if (iy > 0) // bottom: cell is c2
            penalty -= JumpPenaltyCouplingX[cell].transpose() * x.row(cell - nx).transpose();
        if (ix + 1 < nx) // right: cell is c1
            penalty -= JumpPenaltyCouplingY[leftEdge + 1] * x.row(cell + 1).transpose();
        if (iy + 1 < ny) // top: cell is c1
            penalty -= JumpPenaltyCouplingX[cell + nx] * x.row(cell + nx).transpose();
        if (ix > 0) // left: cell is c2
            penalty -= JumpPenaltyCouplingY[leftEdge].transpose() * x.row(cell - 1).transpose();

        y.row(cell) = penaltyCoefficient * penalty.transpose();
        if (massCoefficient != 0.0)
            y.row(cell) += massCoefficient
                * (WeightedDGMassMatrix[cell] * x.row(cell).transpose()).transpose();
    }
}

// LDG: velocity, stress and strain in the same space
template class ParametricDGMomentumMap<3, 3, 3>;
template class ParametricDGMomentumMap<4, 4, 4>;
template class ParametricDGMomentumMap<6, 6, 6>;
template class ParametricDGMomentumMap<9, 9, 9>;
template void ParametricDGMomentumMap<3, 3, 3>::InitializeWeightedDGMassMatrix<1>(
    const DGVector<1>&);
template void ParametricDGMomentumMap<3, 3, 3>::InitializeWeightedDGMassMatrix<3>(
    const DGVector<3>&);
template void ParametricDGMomentumMap<3, 3, 3>::InitializeWeightedDGMassMatrix<4>(
    const DGVector<4>&);
template void ParametricDGMomentumMap<3, 3, 3>::InitializeWeightedDGMassMatrix<6>(
    const DGVector<6>&);
template void ParametricDGMomentumMap<4, 4, 4>::InitializeWeightedDGMassMatrix<1>(
    const DGVector<1>&);
template void ParametricDGMomentumMap<4, 4, 4>::InitializeWeightedDGMassMatrix<3>(
    const DGVector<3>&);
template void ParametricDGMomentumMap<4, 4, 4>::InitializeWeightedDGMassMatrix<4>(
    const DGVector<4>&);
template void ParametricDGMomentumMap<4, 4, 4>::InitializeWeightedDGMassMatrix<6>(
    const DGVector<6>&);
template void ParametricDGMomentumMap<6, 6, 6>::InitializeWeightedDGMassMatrix<1>(
    const DGVector<1>&);
template void ParametricDGMomentumMap<6, 6, 6>::InitializeWeightedDGMassMatrix<3>(
    const DGVector<3>&);
template void ParametricDGMomentumMap<6, 6, 6>::InitializeWeightedDGMassMatrix<4>(
    const DGVector<4>&);
template void ParametricDGMomentumMap<6, 6, 6>::InitializeWeightedDGMassMatrix<6>(
    const DGVector<6>&);
template void ParametricDGMomentumMap<9, 9, 9>::InitializeWeightedDGMassMatrix<1>(
    const DGVector<1>&);
template void ParametricDGMomentumMap<9, 9, 9>::InitializeWeightedDGMassMatrix<3>(
    const DGVector<3>&);
template void ParametricDGMomentumMap<9, 9, 9>::InitializeWeightedDGMassMatrix<4>(
    const DGVector<4>&);
template void ParametricDGMomentumMap<9, 9, 9>::InitializeWeightedDGMassMatrix<6>(
    const DGVector<6>&);

#define INSTANTIATE_VISCOUS_JUMP_PENALTY(DG, DGadvection)                                          \
    template void ParametricDGMomentumMap<DG, DG, DG>::UpdateViscousJumpPenalty<DGadvection>(      \
        const FloatType, const FloatType, const FloatType, const FloatType, const VPParameters&,   \
        const DGVector<DG>&, const DGVector<DG>&, const DGVector<DG>&,                             \
        const DGVector<DGadvection>&, const DGVector<DGadvection>&);
#define INSTANTIATE_VISCOUS_JUMP_PENALTY_ALL(DG)                                                   \
    INSTANTIATE_VISCOUS_JUMP_PENALTY(DG, 1)                                                        \
    INSTANTIATE_VISCOUS_JUMP_PENALTY(DG, 3)                                                        \
    INSTANTIATE_VISCOUS_JUMP_PENALTY(DG, 4)                                                        \
    INSTANTIATE_VISCOUS_JUMP_PENALTY(DG, 6)

INSTANTIATE_VISCOUS_JUMP_PENALTY_ALL(3)
INSTANTIATE_VISCOUS_JUMP_PENALTY_ALL(4)
INSTANTIATE_VISCOUS_JUMP_PENALTY_ALL(6)
INSTANTIATE_VISCOUS_JUMP_PENALTY_ALL(9)

#undef INSTANTIATE_VISCOUS_JUMP_PENALTY_ALL
#undef INSTANTIATE_VISCOUS_JUMP_PENALTY

}
