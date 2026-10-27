/*!
 * @file LDGMomentum.cpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 * Local discontinuous Galerkin (LDG) discretization of the momentum equation with the mEVP
 * solver.
 */

#include "LDGMomentum.hpp"

#include "Interpolations.hpp"
#include "MEVPStressUpdateStep.hpp"
#include "ParametricTools.hpp"
#include "cgVector.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgEdgeTools.hpp"
#include "dgLimit.hpp"
#include "dgPCG.hpp"
#include "include/constants.hpp"

#include <algorithm>
#include <iostream>

namespace Nextsim {

//////////////////////////////////////////////////
// LDG flux terms
//////////////////////////////////////////////////
//
// Interior edge between c1 ("-" side) and c2 ("+" side), normal n from c1 into c2:
//   sigma_hat.n = (0.5-a) sigma1.n + (0.5+a) sigma2.n
//   u_hat       = (0.5+a) u1       + (0.5-a) u2
// Edge integrals: int_e f ds = length * sum_g w_g f(xi_g) (straight edges), PSIe_w<DG,NGP,side>
// holds w_g * phi(Gauss point g of that side).

namespace {

    //! Length and unit normal of a straight edge. The normal is the tangent (end - start)
    //! rotated by +90 degrees (normalSense = +1) or -90 degrees (normalSense = -1).
    struct EdgeGeometry {
        FloatType length;
        FloatType normal_x;
        FloatType normal_y;
    };

    EdgeGeometry edgeGeometry(
        const Eigen::Matrix<FloatType, 2, 2>& ends, const FloatType normalSense)
    {
        const Eigen::Matrix<FloatType, 1, 2> tangent = ends.row(1) - ends.row(0);
        const FloatType length = tangent.norm();
        return { length, -normalSense * tangent(1) / length, normalSense * tangent(0) / length };
    }

    //! Outer normal of the edge on side Side (0 bottom, 1 right, 2 top, 3 left) of cell
    template <int Side> EdgeGeometry edgeOfCell(const ParametricMesh& smesh, const size_t cell)
    {
        const size_t nx = smesh.nx;
        const size_t leftEdge = (cell / nx) * (nx + 1) + cell % nx;
        if constexpr (Side == 0)
            return edgeGeometry(smesh.coordinatesOfEdgeX(cell), -1.0);
        else if constexpr (Side == 1)
            return edgeGeometry(smesh.coordinatesOfEdgeY(leftEdge + 1), -1.0);
        else if constexpr (Side == 2)
            return edgeGeometry(smesh.coordinatesOfEdgeX(cell + nx), +1.0);
        else
            return edgeGeometry(smesh.coordinatesOfEdgeY(leftEdge), +1.0);
    }

    //! Values of a DG field in the NGP Gauss points of side Side of a cell
    template <int DGfield, int NGP, int Side>
    LocalEdgeVector<NGP> traceOnSide(const DGVector<DGfield>& field, const size_t cell)
    {
        if constexpr (Side == 0)
            return bottomedgeofcell<DGfield>(field, cell) * PSIe<EDGEDOFS(DGfield), NGP>;
        else if constexpr (Side == 1)
            return rightedgeofcell<DGfield>(field, cell) * PSIe<EDGEDOFS(DGfield), NGP>;
        else if constexpr (Side == 2)
            return topedgeofcell<DGfield>(field, cell) * PSIe<EDGEDOFS(DGfield), NGP>;
        else
            return leftedgeofcell<DGfield>(field, cell) * PSIe<EDGEDOFS(DGfield), NGP>;
    }

    /*!
     * -<sigma_hat.n, [[v]]>_e for one interior edge, seen by c1 on its side Side1 and by c2 on
     * its side Side2. [[v]] = v1 - v2: subtracted in c1, added in c2.
     */
    template <int DG, int DGstress, int NGP, int Side1, int Side2>
    void sigmaFluxInteriorEdge(const FloatType a, const EdgeGeometry& edge,
        const DGVector<DGstress>& S11, const DGVector<DGstress>& S12, const DGVector<DGstress>& S22,
        DGVector<DG>& vx_residual, DGVector<DG>& vy_residual, const size_t c1, const size_t c2)
    {
        const FloatType nx = edge.normal_x;
        const FloatType ny = edge.normal_y;

        // sigma.n = (S11 nx + S12 ny, S12 nx + S22 ny) on both sides
        const LocalEdgeVector<NGP> sigma1_n_x = traceOnSide<DGstress, NGP, Side1>(S11, c1) * nx
            + traceOnSide<DGstress, NGP, Side1>(S12, c1) * ny;
        const LocalEdgeVector<NGP> sigma1_n_y = traceOnSide<DGstress, NGP, Side1>(S12, c1) * nx
            + traceOnSide<DGstress, NGP, Side1>(S22, c1) * ny;
        const LocalEdgeVector<NGP> sigma2_n_x = traceOnSide<DGstress, NGP, Side2>(S11, c2) * nx
            + traceOnSide<DGstress, NGP, Side2>(S12, c2) * ny;
        const LocalEdgeVector<NGP> sigma2_n_y = traceOnSide<DGstress, NGP, Side2>(S12, c2) * nx
            + traceOnSide<DGstress, NGP, Side2>(S22, c2) * ny;

        const LocalEdgeVector<NGP> sigmahat_n_x
            = edge.length * ((0.5_ft - a) * sigma1_n_x + (0.5_ft + a) * sigma2_n_x);
        const LocalEdgeVector<NGP> sigmahat_n_y
            = edge.length * ((0.5_ft - a) * sigma1_n_y + (0.5_ft + a) * sigma2_n_y);

        vx_residual.row(c1) -= sigmahat_n_x * PSIe_w<DG, NGP, Side1>;
        vx_residual.row(c2) += sigmahat_n_x * PSIe_w<DG, NGP, Side2>;
        vy_residual.row(c1) -= sigmahat_n_y * PSIe_w<DG, NGP, Side1>;
        vy_residual.row(c2) += sigmahat_n_y * PSIe_w<DG, NGP, Side2>;
    }

    //! -<sigma_hat.n, v>_e on the Dirichlet edges of side Side: sigma_hat.n = sigma.n
    template <int DG, int DGstress, int NGP, int Side>
    void sigmaFluxDirichlet(const ParametricMesh& smesh, const DGVector<DGstress>& S11,
        const DGVector<DGstress>& S12, const DGVector<DGstress>& S22, DGVector<DG>& vx_residual,
        DGVector<DG>& vy_residual)
    {
#pragma omp parallel for
        for (size_t i = 0; i < smesh.dirichlet[Side].size(); ++i) {
            const size_t cell = smesh.dirichlet[Side][i];
            if (smesh.landmask[cell] == 0)
                continue;
            const EdgeGeometry edge = edgeOfCell<Side>(smesh, cell);
            const LocalEdgeVector<NGP> sigmahat_n_x = edge.length
                * (traceOnSide<DGstress, NGP, Side>(S11, cell) * edge.normal_x
                    + traceOnSide<DGstress, NGP, Side>(S12, cell) * edge.normal_y);
            const LocalEdgeVector<NGP> sigmahat_n_y = edge.length
                * (traceOnSide<DGstress, NGP, Side>(S12, cell) * edge.normal_x
                    + traceOnSide<DGstress, NGP, Side>(S22, cell) * edge.normal_y);
            vx_residual.row(cell) -= sigmahat_n_x * PSIe_w<DG, NGP, Side>;
            vy_residual.row(cell) -= sigmahat_n_y * PSIe_w<DG, NGP, Side>;
        }
    }

    /*!
     * Adds <u_hat, tau.n>_e of one side of a cell (outer normal n) to E11/E12/E22, multiplied
     * with the local inverse strain mass matrix. tau.n = (tau11 nx + tau12 ny, tau12 nx +
     * tau22 ny): E11 sees u_hat_x nx, E22 sees u_hat_y ny, E12 sees 0.5 (u_hat_x ny + u_hat_y
     * nx). Together with the cell term this gives 0.5 (grad u + grad u^T).
     */
    template <int DG, int DGstress, int DGstrain, int NGP, int Side>
    void addVelocityFlux(const ParametricDGMomentumMap<DG, DGstress, DGstrain>& pmap,
        const EdgeGeometry& edge, const LocalEdgeVector<NGP>& uhat_x,
        const LocalEdgeVector<NGP>& uhat_y, DGVector<DGstrain>& E11, DGVector<DGstrain>& E12,
        DGVector<DGstrain>& E22, const size_t cell)
    {
        const LocalEdgeVector<NGP> uhat_n_11 = edge.length * edge.normal_x * uhat_x;
        const LocalEdgeVector<NGP> uhat_n_22 = edge.length * edge.normal_y * uhat_y;
        const LocalEdgeVector<NGP> uhat_n_12
            = 0.5_ft * edge.length * (edge.normal_y * uhat_x + edge.normal_x * uhat_y);

        E11.row(cell)
            += (uhat_n_11 * PSIe_w<DGstrain, NGP, Side>)*pmap.InverseStrainMassMatrix[cell];
        E22.row(cell)
            += (uhat_n_22 * PSIe_w<DGstrain, NGP, Side>)*pmap.InverseStrainMassMatrix[cell];
        E12.row(cell)
            += (uhat_n_12 * PSIe_w<DGstrain, NGP, Side>)*pmap.InverseStrainMassMatrix[cell];
    }

    //! <u_hat, tau.n>_e for one interior edge, c1 on its side Side1 (normal into c2)
    template <int DG, int DGstress, int DGstrain, int NGP, int Side1, int Side2>
    void velocityFluxInteriorEdge(const ParametricDGMomentumMap<DG, DGstress, DGstrain>& pmap,
        const FloatType a, const EdgeGeometry& edge, const DGVector<DG>& vx, const DGVector<DG>& vy,
        DGVector<DGstrain>& E11, DGVector<DGstrain>& E12, DGVector<DGstrain>& E22, const size_t c1,
        const size_t c2)
    {
        const LocalEdgeVector<NGP> uhat_x = (0.5_ft + a) * traceOnSide<DG, NGP, Side1>(vx, c1)
            + (0.5_ft - a) * traceOnSide<DG, NGP, Side2>(vx, c2);
        const LocalEdgeVector<NGP> uhat_y = (0.5_ft + a) * traceOnSide<DG, NGP, Side1>(vy, c1)
            + (0.5_ft - a) * traceOnSide<DG, NGP, Side2>(vy, c2);

        // the normal points out of c1 and into c2
        const EdgeGeometry inverted = { edge.length, -edge.normal_x, -edge.normal_y };
        addVelocityFlux<DG, DGstress, DGstrain, NGP, Side1>(
            pmap, edge, uhat_x, uhat_y, E11, E12, E22, c1);
        addVelocityFlux<DG, DGstress, DGstrain, NGP, Side2>(
            pmap, inverted, uhat_x, uhat_y, E11, E12, E22, c2);
    }

    //! <u, tau.n>_e on the free edges of side Side (u_hat = trace of the cell)
    template <int DG, int DGstress, int DGstrain, int NGP, int Side>
    void velocityFluxFree(const ParametricMesh& smesh,
        const ParametricDGMomentumMap<DG, DGstress, DGstrain>& pmap,
        const std::vector<size_t>& cells, const DGVector<DG>& vx, const DGVector<DG>& vy,
        DGVector<DGstrain>& E11, DGVector<DGstrain>& E12, DGVector<DGstrain>& E22)
    {
#pragma omp parallel for
        for (size_t i = 0; i < cells.size(); ++i) {
            const size_t cell = cells[i];
            addVelocityFlux<DG, DGstress, DGstrain, NGP, Side>(pmap, edgeOfCell<Side>(smesh, cell),
                traceOnSide<DG, NGP, Side>(vx, cell), traceOnSide<DG, NGP, Side>(vy, cell), E11,
                E12, E22, cell);
        }
    }

} // anonymous namespace

std::array<std::vector<size_t>, 4> freeEdgesOfMesh(const ParametricMesh& smesh)
{
    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;
    std::array<std::vector<size_t>, 4> freeEdges;
    auto add = [&](int side, size_t cell) {
        if (smesh.landmask[cell] == 0)
            return;
        const auto& d = smesh.dirichlet[side];
        if (std::find(d.begin(), d.end(), cell) == d.end())
            freeEdges[side].push_back(cell);
    };
    for (size_t ix = 0; ix < nx; ++ix) {
        add(0, ix);
        add(2, (ny - 1) * nx + ix);
    }
    for (size_t iy = 0; iy < ny; ++iy) {
        add(1, iy * nx + nx - 1);
        add(3, iy * nx);
    }
    return freeEdges;
}

template <int DG, int DGstress>
void ldgSigmaFlux(const ParametricMesh& smesh, const FloatType a, const DGVector<DGstress>& S11,
    const DGVector<DGstress>& S12, const DGVector<DGstress>& S22, DGVector<DG>& vx_residual,
    DGVector<DG>& vy_residual)
{
    // the integrand sigma * v has degree DGDEGREE(DG)+DGDEGREE(DGstress) along the edge
    constexpr int NGP = NGAUSSFORDEGREE1D(DGDEGREE(DG) + DGDEGREE(DGstress));
    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;

    // Interior X-edges, one column of cells per thread: no two threads write the same cell
#pragma omp parallel for
    for (size_t ix = 0; ix < nx; ++ix)
        for (size_t iy = 0; iy + 1 < ny; ++iy) {
            const size_t cellBelow = iy * nx + ix;
            const size_t cellAbove = cellBelow + nx;
            if (smesh.landmask[cellBelow] == 0 || smesh.landmask[cellAbove] == 0)
                continue;
            // bottom edge of cellAbove, normal from cellBelow into cellAbove
            const EdgeGeometry edge = edgeGeometry(smesh.coordinatesOfEdgeX(cellAbove), +1.0);
            sigmaFluxInteriorEdge<DG, DGstress, NGP, 2, 0>(
                a, edge, S11, S12, S22, vx_residual, vy_residual, cellBelow, cellAbove);
        }

    // Interior Y-edges, one row of cells per thread
#pragma omp parallel for
    for (size_t iy = 0; iy < ny; ++iy)
        for (size_t ix = 0; ix + 1 < nx; ++ix) {
            const size_t cellLeft = iy * nx + ix;
            const size_t cellRight = cellLeft + 1;
            if (smesh.landmask[cellLeft] == 0 || smesh.landmask[cellRight] == 0)
                continue;
            // right edge of cellLeft, normal from cellLeft into cellRight
            const EdgeGeometry edge
                = edgeGeometry(smesh.coordinatesOfEdgeY(iy * (nx + 1) + ix + 1), -1.0);
            sigmaFluxInteriorEdge<DG, DGstress, NGP, 1, 3>(
                a, edge, S11, S12, S22, vx_residual, vy_residual, cellLeft, cellRight);
        }

    // Dirichlet edges (the sides one after the other: a corner cell is in two lists)
    sigmaFluxDirichlet<DG, DGstress, NGP, 0>(smesh, S11, S12, S22, vx_residual, vy_residual);
    sigmaFluxDirichlet<DG, DGstress, NGP, 1>(smesh, S11, S12, S22, vx_residual, vy_residual);
    sigmaFluxDirichlet<DG, DGstress, NGP, 2>(smesh, S11, S12, S22, vx_residual, vy_residual);
    sigmaFluxDirichlet<DG, DGstress, NGP, 3>(smesh, S11, S12, S22, vx_residual, vy_residual);

    // free edges: sigma_hat.n = 0, no contribution
}

template <int DG, int DGstress, int DGstrain>
void ldgVelocityFlux(const ParametricMesh& smesh,
    const ParametricDGMomentumMap<DG, DGstress, DGstrain>& pmap, const FloatType a,
    const std::array<std::vector<size_t>, 4>& freeEdges, const DGVector<DG>& vx,
    const DGVector<DG>& vy, DGVector<DGstrain>& E11, DGVector<DGstrain>& E12,
    DGVector<DGstrain>& E22)
{
    // the integrand u * tau has degree DGDEGREE(DG)+DGDEGREE(DGstrain) along the edge
    constexpr int NGP = NGAUSSFORDEGREE1D(DGDEGREE(DG) + DGDEGREE(DGstrain));
    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;

    // Interior X-edges, one column of cells per thread
#pragma omp parallel for
    for (size_t ix = 0; ix < nx; ++ix)
        for (size_t iy = 0; iy + 1 < ny; ++iy) {
            const size_t cellBelow = iy * nx + ix;
            const size_t cellAbove = cellBelow + nx;
            if (smesh.landmask[cellBelow] == 0 || smesh.landmask[cellAbove] == 0)
                continue;
            // bottom edge of cellAbove, normal from cellBelow into cellAbove
            const EdgeGeometry edge = edgeGeometry(smesh.coordinatesOfEdgeX(cellAbove), +1.0);
            velocityFluxInteriorEdge<DG, DGstress, DGstrain, NGP, 2, 0>(
                pmap, a, edge, vx, vy, E11, E12, E22, cellBelow, cellAbove);
        }

    // Interior Y-edges, one row of cells per thread
#pragma omp parallel for
    for (size_t iy = 0; iy < ny; ++iy)
        for (size_t ix = 0; ix + 1 < nx; ++ix) {
            const size_t cellLeft = iy * nx + ix;
            const size_t cellRight = cellLeft + 1;
            if (smesh.landmask[cellLeft] == 0 || smesh.landmask[cellRight] == 0)
                continue;
            // right edge of cellLeft, normal from cellLeft into cellRight
            const EdgeGeometry edge
                = edgeGeometry(smesh.coordinatesOfEdgeY(iy * (nx + 1) + ix + 1), -1.0);
            velocityFluxInteriorEdge<DG, DGstress, DGstrain, NGP, 1, 3>(
                pmap, a, edge, vx, vy, E11, E12, E22, cellLeft, cellRight);
        }

    // Dirichlet edges: u_hat = 0, no contribution. Free edges: u_hat = u.
    velocityFluxFree<DG, DGstress, DGstrain, NGP, 0>(
        smesh, pmap, freeEdges[0], vx, vy, E11, E12, E22);
    velocityFluxFree<DG, DGstress, DGstrain, NGP, 1>(
        smesh, pmap, freeEdges[1], vx, vy, E11, E12, E22);
    velocityFluxFree<DG, DGstress, DGstrain, NGP, 2>(
        smesh, pmap, freeEdges[2], vx, vy, E11, E12, E22);
    velocityFluxFree<DG, DGstress, DGstrain, NGP, 3>(
        smesh, pmap, freeEdges[3], vx, vy, E11, E12, E22);
}

//////////////////////////////////////////////////
// LDGMomentum
//////////////////////////////////////////////////

template <int DG, int DGadvection>
LDGMomentum<DG, DGadvection>::LDGMomentum(
    const ParametricMesh& smesh_, const FloatType a_, const FloatType b)
    : pmap(smesh_)
    , smesh(smesh_)
    , a(a_)
    , freeEdges(freeEdgesOfMesh(smesh_))
{
    if (!(smesh.nelements > 0)) {
        std::cerr << "LDGMomentum: The mesh has to be initialized first!" << std::endl;
        abort();
    }
    if (smesh.CoordinateSystem != CARTESIAN) {
        std::cerr << "LDGMomentum: only CARTESIAN coordinates are supported" << std::endl;
        abort();
    }
    if (!smesh.periodic.empty()) {
        std::cerr << "LDGMomentum: periodic boundaries are not supported" << std::endl;
        abort();
    }

    for (DGVector<DG>* v :
        { &vx, &vy, &vx0, &vy0, &ox, &oy, &ax, &ay, &E11, &E12, &E22, &S11, &S12, &S22 }) {
        v->resize_by_mesh(smesh);
        v->setZero();
    }
    H.resize_by_mesh(smesh);
    A.resize_by_mesh(smesh);
    H.setZero();
    A.setZero();
    sshGradX.setZero(smesh.nelements, GAUSSPOINTS(DG));
    sshGradY.setZero(smesh.nelements, GAUSSPOINTS(DG));

    pmap.InitializeDGDivMatrices();
    pmap.InitializeJwPSI();
    pmap.InitializeIMJwPSI();
    pmap.InitializeJumpPenalty(a, b);
}

template <int DG, int DGadvection>
void LDGMomentum<DG, DGadvection>::setSeaSurfaceHeight(const DGVector<1>& ssh)
{
    // continuous (CG1) representation of the sea surface height
    CGVector<1> cgssh(smesh);
    Interpolations::DG2CG(smesh, cgssh, ssh);

    constexpr int Q = GAUSSPOINTS1D(DG);
    const size_t nx = smesh.nx;
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        const size_t n0 = (eid / nx) * (nx + 1) + eid % nx; // lower left node
        const Eigen::Matrix<FloatType, 1, 4> local
            = { cgssh(n0), cgssh(n0 + 1), cgssh(n0 + nx + 1), cgssh(n0 + nx + 2) };
        // derivatives in the reference coordinates
        const Eigen::Matrix<FloatType, 1, Q * Q> fx = local * PHIx<1, Q>;
        const Eigen::Matrix<FloatType, 1, Q * Q> fy = local * PHIy<1, Q>;
        const Eigen::Matrix<FloatType, 2, Q * Q> dxT = ParametricTools::dxT<Q>(smesh, eid);
        const Eigen::Matrix<FloatType, 2, Q * Q> dyT = ParametricTools::dyT<Q>(smesh, eid);
        const Eigen::Matrix<FloatType, 1, Q * Q> J = ParametricTools::J<Q>(smesh, eid);
        // grad f = dT^{-T} grad_ref f = 1/J [ dyT2 fx - dxT2 fy, -dyT1 fx + dxT1 fy ]
        sshGradX.row(eid)
            = ((dyT.row(1).array() * fx.array() - dxT.row(1).array() * fy.array()) / J.array())
                  .matrix();
        sshGradY.row(eid)
            = ((dxT.row(0).array() * fy.array() - dyT.row(0).array() * fx.array()) / J.array())
                  .matrix();
    }
}

template <int DG, int DGadvection>
void LDGMomentum<DG, DGadvection>::prepareIteration(const VPParameters& params,
    const DGVector<DGadvection>& Hin, const DGVector<DGadvection>& Ain, const FloatType beta,
    const FloatType dt)
{
    vx0 = vx;
    vy0 = vy;

    // avoid H = 0 for the invertibility of the mass matrix
    H = Hin;
    LimitMin(H, 1.e-4);
    A = Ain;

    // also refreshes pmap.InverseWeightedDGMassMatrix, used as the PCG preconditioner
    pmap.template InitializeWeightedDGMassMatrix<DGadvection>(H);

    velocityMassFactor = (1.0_ft + beta) * params.rhoIce;
    velocityDt = dt;
}

template <int DG, int DGadvection> void LDGMomentum<DG, DGadvection>::projectVelocityToStrain()
{
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        if (smesh.landmask[eid] == 0) {
            E11.row(eid).setZero();
            E12.row(eid).setZero();
            E22.row(eid).setZero();
            continue;
        }
        E11.row(eid) = -(pmap.iMgradTauX[eid] * vx.row(eid).transpose());
        E22.row(eid) = -(pmap.iMgradTauY[eid] * vy.row(eid).transpose());
        E12.row(eid) = -0.5_ft
            * (pmap.iMgradTauY[eid] * vx.row(eid).transpose()
                + pmap.iMgradTauX[eid] * vy.row(eid).transpose());
    }
    ldgVelocityFlux<DG, DG, DG>(smesh, pmap, a, freeEdges, vx, vy, E11, E12, E22);
}

template <int DG, int DGadvection>
void LDGMomentum<DG, DGadvection>::mEVPStep(
    const VPParameters& params, const FloatType alpha, const FloatType beta, const FloatType dt)
{
    constexpr int NGP1 = GAUSSPOINTS1D(DG);
    constexpr int NGP = GAUSSPOINTS(DG);

    // (a) stress update
    mEVPStressUpdate<DGadvection, DG>(
        params, smesh, pmap.iMJwPSI, S11, S12, S22, E11, E12, E22, H, A, alpha);

    // (b) velocity: ((1+beta) rho_ice M_H + dt J) u = F
    const FloatType FOcean = params.COcean * params.rhoOcean;
    const FloatType FAtm = params.CAtm * params.rhoAtm;

    DGVector<DG> Fx(smesh), Fy(smesh);
    Fx.setZero();
    Fy.setZero();

    // forcing, relaxation and cell term of (sigma, grad v)
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        if (smesh.landmask[eid] == 0)
            continue;

        const Eigen::Matrix<FloatType, 1, NGP> vx_g = vx.row(eid) * PSI<DG, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> vy_g = vy.row(eid) * PSI<DG, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> ox_g = ox.row(eid) * PSI<DG, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> oy_g = oy.row(eid) * PSI<DG, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> ax_g = ax.row(eid) * PSI<DG, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> ay_g = ay.row(eid) * PSI<DG, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> H_g = H.row(eid) * PSI<DGadvection, NGP1>;
        const Eigen::Matrix<FloatType, 1, NGP> A_g
            = (A.row(eid) * PSI<DGadvection, NGP1>).array().max(0.0_ft).min(1.0_ft).matrix();

        const Eigen::Matrix<FloatType, 1, NGP> absatm
            = (ax_g.array().square() + ay_g.array().square()).sqrt().matrix();
        const Eigen::Matrix<FloatType, 1, NGP> absocn
            = ((vx_g - ox_g).array().square() + (vy_g - oy_g).array().square()).sqrt().matrix();

        // wind and ocean drag, Coriolis, sea surface tilt (explicit)
        const Eigen::Matrix<FloatType, 1, NGP> Force_x = (A_g.array()
                * (FAtm * absatm.array() * ax_g.array()
                    + FOcean * absocn.array() * (ox_g.array() - vx_g.array()))
            + params.rhoIce * H_g.array()
                * (params.fc * vy_g.array() - PhysicalConstants::g * sshGradX.row(eid).array()))
                                                             .matrix();
        const Eigen::Matrix<FloatType, 1, NGP> Force_y = (A_g.array()
                * (FAtm * absatm.array() * ay_g.array()
                    + FOcean * absocn.array() * (oy_g.array() - vy_g.array()))
            - params.rhoIce * H_g.array()
                * (params.fc * vx_g.array() + PhysicalConstants::g * sshGradY.row(eid).array()))
                                                             .matrix();

        const Eigen::Matrix<FloatType, DG, 1> Load_x = pmap.JwPSI[eid] * Force_x.transpose();
        const Eigen::Matrix<FloatType, DG, 1> Load_y = pmap.JwPSI[eid] * Force_y.transpose();

        // relaxation rho_ice M_H (beta u^k + u^n)
        const Eigen::Matrix<FloatType, DG, 1> Relax_x = params.rhoIce
            * pmap.WeightedDGMassMatrix[eid] * (beta * vx.row(eid) + vx0.row(eid)).transpose();
        const Eigen::Matrix<FloatType, DG, 1> Relax_y = params.rhoIce
            * pmap.WeightedDGMassMatrix[eid] * (beta * vy.row(eid) + vy0.row(eid)).transpose();

        // cell term of (sigma, grad v)
        const Eigen::Matrix<FloatType, DG, 1> SigmaCell_x
            = pmap.dgDivS1[eid] * S11.row(eid).transpose()
            + pmap.dgDivS2[eid] * S12.row(eid).transpose();
        const Eigen::Matrix<FloatType, DG, 1> SigmaCell_y
            = pmap.dgDivS1[eid] * S12.row(eid).transpose()
            + pmap.dgDivS2[eid] * S22.row(eid).transpose();

        Fx.row(eid) = (dt * Load_x + Relax_x - dt * SigmaCell_x).transpose();
        Fy.row(eid) = (dt * Load_y + Relax_y - dt * SigmaCell_y).transpose();
    }

    // edge term -<sigma_hat.n, [[v]]> of the known stress
    DGVector<DG> FluxSigma_x(smesh), FluxSigma_y(smesh);
    FluxSigma_x.setZero();
    FluxSigma_y.setZero();
    ldgSigmaFlux<DG, DG>(smesh, a, S11, S12, S22, FluxSigma_x, FluxSigma_y);
    Fx -= dt * FluxSigma_x;
    Fy -= dt * FluxSigma_y;

    // matrix-free PCG, preconditioned with the cell mass blocks. vx/vy are the initial guess.
    auto applyA = [&](const DGVector<DG>& x, DGVector<DG>& y) {
        pmap.ApplyMassPlusJumpPenalty(velocityMassFactor, velocityDt, x, y);
    };
    auto applyPrec = [&](const DGVector<DG>& r, DGVector<DG>& z) {
#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid)
            z.row(eid)
                = (pmap.InverseWeightedDGMassMatrix[eid] * r.row(eid).transpose()).transpose()
                / velocityMassFactor;
    };
    FloatType err_x, err_y;
    const int it_x = PCG<DG>(applyA, applyPrec, Fx, vx, pcgTolerance, pcgMaxIterations, err_x);
    const int it_y = PCG<DG>(applyA, applyPrec, Fy, vy, pcgTolerance, pcgMaxIterations, err_y);
    pcgIterations = std::max(it_x, it_y);
    if (it_x >= pcgMaxIterations || it_y >= pcgMaxIterations)
        std::cerr << "LDGMomentum: PCG not converged (vx err " << err_x << ", vy err " << err_y
                  << ")" << std::endl;

    // velocity is zero on land
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        if (smesh.landmask[eid] == 0) {
            vx.row(eid).setZero();
            vy.row(eid).setZero();
        }

    // (c) strain
    projectVelocityToStrain();
}

// velocity, stress and strain in DG3, DG4, DG6, DG9; advection in DG1, DG3, DG4, DG6
#define INSTANTIATE_LDG(DG)                                                                        \
    template class LDGMomentum<DG, 1>;                                                             \
    template class LDGMomentum<DG, 3>;                                                             \
    template class LDGMomentum<DG, 4>;                                                             \
    template class LDGMomentum<DG, 6>;                                                             \
    template void ldgSigmaFlux<DG, DG>(const ParametricMesh&, const FloatType,                     \
        const DGVector<DG>&, const DGVector<DG>&, const DGVector<DG>&, DGVector<DG>&,              \
        DGVector<DG>&);                                                                            \
    template void ldgVelocityFlux<DG, DG, DG>(const ParametricMesh&,                               \
        const ParametricDGMomentumMap<DG, DG, DG>&, const FloatType,                               \
        const std::array<std::vector<size_t>, 4>&, const DGVector<DG>&, const DGVector<DG>&,       \
        DGVector<DG>&, DGVector<DG>&, DGVector<DG>&);

INSTANTIATE_LDG(3)
INSTANTIATE_LDG(4)
INSTANTIATE_LDG(6)
INSTANTIATE_LDG(9)

#undef INSTANTIATE_LDG

} /* namespace Nextsim */
