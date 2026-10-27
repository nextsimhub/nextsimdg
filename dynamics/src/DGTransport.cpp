/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#include "DGTransport.hpp"
#include "Interpolations.hpp"
#include "ParametricTools.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgEdgeTools.hpp"

namespace Nextsim {

// leftedgeofcell, rightedgeofcell, bottomedgeofcell, topedgeofcell: see dgEdgeTools.hpp

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Computes the normal velocity on each edge
// The velocity is scaled with the length of the edge, therefore
// the length-element J is already included!
//
// The way to compute the normal vector does not differ between
// the Cartesian and the Spherical version. This is, since the normal
// is not the actual normal but a scaled version of it. The weight
// that is missing is corrected in the two functions
// edge_term_X and edge_term_Y. Here, we omit the scaling with
// cos(theta) and simple use ds = |J| as length scale
// (instead of ds = (d_theta^2 + cos(theta)^2 d_phi^2)^1/2
//
template <int DG> void DGTransport<DG>::reinitnormalvelocity()
{
    // average the velocity to the Y-edges
    normalvel_Y.zero(); // < Parallelize
    normalvel_X.zero();

#pragma omp parallel for
    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        //   |     |
        // --*-----*--
        //  ey  cy |
        //   |     |
        // -ey-----*--
        //   |     |
        size_t ey = iy * (smesh.nx + 1); // first edge-index and node-index in row
        size_t cy = iy * smesh.nx; // first cell index in row

        for (size_t ix = 0; ix < smesh.nx; ++ix, ++ey, ++cy) {
            if (smesh.landmask[cy] == 0) // skip all land elements
                continue;

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<FloatType, 1, 2> tangent_left
                = smesh.edgevector(ey, ey + smesh.nx + 1);
            normalvel_Y.row(ey) += 0.5
                * (tangent_left(0, 1) * leftedgeofcell<DG>(velx, cy)
                    - tangent_left(0, 0) * leftedgeofcell<DG>(vely, cy));

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<FloatType, 1, 2> tangent_right
                = smesh.edgevector(ey + 1, ey + smesh.nx + 2);
            normalvel_Y.row(ey + 1) += 0.5
                * (tangent_right(0, 1) * rightedgeofcell<DG>(velx, cy)
                    - tangent_right(0, 0) * rightedgeofcell<DG>(vely, cy));
        }
        // we need an adjustment along the boundaries.. This is done later on.
    }

#pragma omp parallel for
    for (size_t ix = 0; ix < smesh.nx; ++ix) {
        //   |     |
        // --*-----*--
        //   |  cx |
        //   |     |
        // -nx-ex--*--
        //   |     |

        size_t cx = ix; // first edge-index and cell-index
        size_t nx = ix; // first cell index in row

        for (size_t iy = 0; iy < smesh.ny; ++iy, cx += smesh.nx, nx += smesh.nx + 1) {
            if (smesh.landmask[cx] == 0) // skip land elements
                continue;

            // un-normed tangent vector of bottom edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<FloatType, 1, 2> tangent_bottom = smesh.edgevector(nx, nx + 1);

            normalvel_X.row(cx) += 0.5
                * (-tangent_bottom(0, 1) * bottomedgeofcell<DG>(velx, cx)
                    + tangent_bottom(0, 0) * bottomedgeofcell<DG>(vely, cx));

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<FloatType, 1, 2> tangent_top
                = smesh.edgevector(nx + smesh.nx + 1, nx + smesh.nx + 2);

            normalvel_X.row(cx + smesh.nx) += 0.5
                * (-tangent_top(0, 1) * topedgeofcell<DG>(velx, cx)
                    + tangent_top(0, 0) * topedgeofcell<DG>(vely, cx));
        }
    }

    // TR 07.04.2025
    // correction of the normal velocity is needed on the mesh-boundary, where each edge
    // only has one neighbor. Above, we sum both sides, weighted with 0.5. Therefore,
    // we must scale the outer edges with 2. This is only used for inflow/outflow
    // Neumann boundaries.
#pragma omp parallel for
    for (size_t i = 0; i < smesh.nx; ++i) {
        normalvel_X.row(i) *= 2.0; // bottom
        normalvel_X.row(i + smesh.nx * smesh.ny) *= 2.0; // top
    }
#pragma omp parallel for
    for (size_t i = 0; i < smesh.ny; ++i) {
        normalvel_Y.row(i * (smesh.nx + 1)) *= 2.0; // left
        normalvel_Y.row(i * (smesh.nx + 1) + smesh.nx) *= 2.0; // right
    }
}

////////////////////////////////////////////////// PREPARE

/*!
 * Prepares the advection step:
 * - interpolates CG velocity to DG
 * - initializes normal velocity on the edges
 */
template <int DG>
template <int CG>
void DGTransport<DG>::prepareAdvection(const CGVector<CG>& cg_vx, const CGVector<CG>& cg_vy)
{
    Nextsim::Interpolations::CG2DG(smesh, GetVx(), cg_vx);
    Nextsim::Interpolations::CG2DG(smesh, GetVy(), cg_vy);
    reinitnormalvelocity();
}

/*!
 * Prepares the advection step with a DG velocity (LDG momentum), see DGTransport.hpp.
 * velx/vely hold the velocity projected onto DG (the basis is hierarchical and orthogonal),
 * they are not used by the transport with the full velocity.
 */
template <int DG>
template <int DGm>
void DGTransport<DG>::prepareAdvectionDG(
    const DGVector<DGm>& dg_vx, const DGVector<DGm>& dg_vy, const FloatType a)
{
    if (!smesh.periodic.empty()) {
        std::cerr << "DGTransport::prepareAdvectionDG: periodic boundaries not supported"
                  << std::endl;
        abort();
    }
    if (cellTermFullX.empty())
        InitializeCellTermsFull();

    constexpr int minDG = std::min(DG, DGm);
    velx.setZero();
    vely.setZero();
    velx.leftCols(minDG) = dg_vx.leftCols(minDG);
    vely.leftCols(minDG) = dg_vy.leftCols(minDG);

    // cell term: full velocity in the NGP_CELL^2 Gauss points of each cell
    cellVel_X.resize(smesh.nelements, NGP_CELL * NGP_CELL);
    cellVel_Y.resize(smesh.nelements, NGP_CELL * NGP_CELL);
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        cellVel_X.row(eid) = dg_vx.row(eid) * PSI<DGm, NGP_CELL>;
        cellVel_Y.row(eid) = dg_vy.row(eid) * PSI<DGm, NGP_CELL>;
    }

    constexpr int NGP = NGP_EDGE;
    const size_t nx = smesh.nx;
    const size_t ny = smesh.ny;

    // 0 on edges touching land
    edgeNormalVel_X.setZero(nx * (ny + 1), NGP);
    edgeNormalVel_Y.setZero((nx + 1) * ny, NGP);

    // trace of the velocity on side 'side' of 'cell' in the Gauss points
    auto trace = [&](const DGVector<DGm>& v, size_t cell, int side) -> LocalEdgeVector<NGP> {
        if (side == 0)
            return bottomedgeofcell<DGm>(v, cell) * PSIe<EDGEDOFS(DGm), NGP>;
        else if (side == 1)
            return rightedgeofcell<DGm>(v, cell) * PSIe<EDGEDOFS(DGm), NGP>;
        else if (side == 2)
            return topedgeofcell<DGm>(v, cell) * PSIe<EDGEDOFS(DGm), NGP>;
        else
            return leftedgeofcell<DGm>(v, cell) * PSIe<EDGEDOFS(DGm), NGP>;
    };

    // X-edges (horizontal), normal (-ty, tx) to +y. Below the edge c1 (top side), above c2
    // (bottom side). Edge iy*nx+ix is the bottom edge of cell (ix,iy).
#pragma omp parallel for
    for (size_t iy = 0; iy <= ny; ++iy)
        for (size_t ix = 0; ix < nx; ++ix) {
            const bool hasBelow = (iy > 0) && smesh.landmask[(iy - 1) * nx + ix];
            const bool hasAbove = (iy < ny) && smesh.landmask[iy * nx + ix];
            const bool landBelow = (iy > 0) && !smesh.landmask[(iy - 1) * nx + ix];
            const bool landAbove = (iy < ny) && !smesh.landmask[iy * nx + ix];
            if (landBelow || landAbove || (!hasBelow && !hasAbove))
                continue;
            LocalEdgeVector<NGP> ux, uy;
            if (hasBelow && hasAbove) {
                ux = (0.5_ft + a) * trace(dg_vx, (iy - 1) * nx + ix, 2)
                    + (0.5_ft - a) * trace(dg_vx, iy * nx + ix, 0);
                uy = (0.5_ft + a) * trace(dg_vy, (iy - 1) * nx + ix, 2)
                    + (0.5_ft - a) * trace(dg_vy, iy * nx + ix, 0);
            } else if (hasBelow) { // top boundary
                ux = trace(dg_vx, (iy - 1) * nx + ix, 2);
                uy = trace(dg_vy, (iy - 1) * nx + ix, 2);
            } else { // bottom boundary
                ux = trace(dg_vx, iy * nx + ix, 0);
                uy = trace(dg_vy, iy * nx + ix, 0);
            }
            const size_t node = iy * (nx + 1) + ix;
            const Eigen::Matrix<FloatType, 1, 2> t = smesh.edgevector(node, node + 1);
            edgeNormalVel_X.row(iy * nx + ix) = -t(0, 1) * ux + t(0, 0) * uy;
        }

    // Y-edges (vertical), normal (ty, -tx) to +x. Left of the edge c1 (right side), right c2
    // (left side). Edge iy*(nx+1)+ix is the left edge of cell (ix,iy).
#pragma omp parallel for
    for (size_t iy = 0; iy < ny; ++iy)
        for (size_t ix = 0; ix <= nx; ++ix) {
            const bool hasLeft = (ix > 0) && smesh.landmask[iy * nx + ix - 1];
            const bool hasRight = (ix < nx) && smesh.landmask[iy * nx + ix];
            const bool landLeft = (ix > 0) && !smesh.landmask[iy * nx + ix - 1];
            const bool landRight = (ix < nx) && !smesh.landmask[iy * nx + ix];
            if (landLeft || landRight || (!hasLeft && !hasRight))
                continue;
            LocalEdgeVector<NGP> ux, uy;
            if (hasLeft && hasRight) {
                ux = (0.5_ft + a) * trace(dg_vx, iy * nx + ix - 1, 1)
                    + (0.5_ft - a) * trace(dg_vx, iy * nx + ix, 3);
                uy = (0.5_ft + a) * trace(dg_vy, iy * nx + ix - 1, 1)
                    + (0.5_ft - a) * trace(dg_vy, iy * nx + ix, 3);
            } else if (hasLeft) { // right boundary
                ux = trace(dg_vx, iy * nx + ix - 1, 1);
                uy = trace(dg_vy, iy * nx + ix - 1, 1);
            } else { // left boundary
                ux = trace(dg_vx, iy * nx + ix, 3);
                uy = trace(dg_vy, iy * nx + ix, 3);
            }
            const size_t ie = iy * (nx + 1) + ix;
            const Eigen::Matrix<FloatType, 1, 2> t = smesh.edgevector(ie, ie + nx + 1);
            edgeNormalVel_Y.row(ie) = t(0, 1) * ux - t(0, 0) * uy;
        }

    projectVelocity = false;
}

////////////////////////////////////////////////// CELL TERM

template <>
void DGTransport<1>::cell_term(const ParametricMesh& smesh, FloatType dt, DGVector<1>& phiup,
    const DGVector<1>& phi, const DGVector<1>& vx, const DGVector<1>& vy, const size_t ic)
{
}

template <int DG>
void DGTransport<DG>::cell_term(const ParametricMesh& smesh, FloatType dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const DGVector<DG>& vx, const DGVector<DG>& vy, const size_t eid)
{
    if (smesh.landmask[eid] == 0) {
        // TODO remove comments that are here to validate correctness of landmask
        // std::cout << "element id " << eid  << " lm " << smesh.landmask[eid] << " " <<
        // phi.row(eid) << std::endl;
        return;
    } else {
        // std::cout << "element id " << eid  << " lm " << smesh.landmask[eid] << " " <<
        // phi.row(eid) << std::endl;
    }
    const Eigen::Matrix<FloatType, 1, GAUSSPOINTS(DG)> vx_gauss
        = vx.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>; //!< velocity in GP
    const Eigen::Matrix<FloatType, 1, GAUSSPOINTS(DG)> vy_gauss
        = vy.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>;
    const Eigen::Matrix<FloatType, 1, GAUSSPOINTS(DG)> phi_gauss
        = (phi.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>).array();

    phiup.row(eid) += dt
        * (parammap.AdvectionCellTermX[eid].array().rowwise() * vx_gauss.array()
            + parammap.AdvectionCellTermY[eid].array().rowwise() * vy_gauss.array())
              .matrix()
        * phi_gauss.transpose();
}
////////////////////////////////////////////////// BOUNDARY HANDLING

template <int DG>
void boundary_lower(const ParametricMesh& smesh, const FloatType dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c,
    const size_t e)
{
    // GP = DGEDGE
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    // block<1, 2>(e, 0) * PSIe<2,2>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp
        = ((bottomedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array()
            * (-vel_gauss).cwiseMax(0).array());
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 0>;
}
template <int DG>
void boundary_upper(const ParametricMesh& smesh, const FloatType dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c,
    const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp
        = ((topedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array()
            * vel_gauss.cwiseMax(0).array());
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 2>;
}
template <int DG>
void boundary_left(const ParametricMesh& smesh, const FloatType dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c,
    const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp
        = ((leftedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array()
            * (-vel_gauss).cwiseMax(0).array());
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 3>;
}
template <int DG>
void boundary_right(const ParametricMesh& smesh, const FloatType dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c,
    const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp
        = ((rightedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array()
            * vel_gauss.cwiseMax(0).array());
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 1>;
}

////////////////////////////////////////////////// EDGE TERMS

template <>
inline void DGTransport<1>::edge_term_X(const ParametricMesh& smesh, const FloatType dt,
    DGVector<1>& phiup, const DGVector<1>& phi, // DG0 (1)
    const EdgeVector<1>& normalvel_X, const size_t c1, const size_t c2, const size_t ie)
{
    if (smesh.landmask[c1] == 0)
        return;
    if (smesh.landmask[c2] == 0)
        return;

    FloatType bottom = phi(c1, 0);
    FloatType top = phi(c2, 0);
    FloatType vel = normalvel_X(ie, 0);

    phiup(c1, 0) -= dt * (std::max(vel, 0.0_ft) * bottom + std::min(vel, 0.0_ft) * top);
    phiup(c2, 0) += dt * (std::max(vel, 0.0_ft) * bottom + std::min(vel, 0.0_ft) * top);
}
template <>
inline void DGTransport<1>::edge_term_Y(const ParametricMesh& smesh, const FloatType dt,
    DGVector<1>& phiup, const DGVector<1>& phi, // DG0 (1)
    const EdgeVector<1>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie)
{
    if (smesh.landmask[c1] == 0)
        return;
    if (smesh.landmask[c2] == 0)
        return;

    FloatType left = phi(c1, 0);
    FloatType right = phi(c2, 0);
    FloatType vel = normalvel_Y(ie, 0);

    phiup(c1, 0) -= dt * (std::max(vel, 0.0_ft) * left + std::min(vel, 0.0_ft) * right);
    phiup(c2, 0) += dt * (std::max(vel, 0.0_ft) * left + std::min(vel, 0.0_ft) * right);
}

template <int DG>
inline void DGTransport<DG>::edge_term_X(const ParametricMesh& smesh, const FloatType dt,
    DGVector<DG>& phiup, const DGVector<DG>& phi, // DG1 (3)
    const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c1, const size_t c2, const size_t ie)
{
    if (smesh.landmask[c1] == 0)
        return;
    if (smesh.landmask[c2] == 0)
        return;

    const LocalEdgeVector<GAUSSPOINTS1D(DG)> vel_gauss
        = normalvel_X.row(ie) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>;

    const LocalEdgeVector<GAUSSPOINTS1D(DG)> tmp = (vel_gauss.cwiseMax(0).array()
            * (topedgeofcell<DG>(phi, c1) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>).array()
        + vel_gauss.cwiseMin(0).array()
            * (bottomedgeofcell<DG>(phi, c2) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>).array());

    phiup.row(c1) -= dt * tmp * PSIe_w<DG, GAUSSPOINTS1D(DG), 2>;
    phiup.row(c2) += dt * tmp * PSIe_w<DG, GAUSSPOINTS1D(DG), 0>;
}

template <int DG>
inline void DGTransport<DG>::edge_term_Y(const ParametricMesh& smesh, const FloatType dt,
    DGVector<DG>& phiup, const DGVector<DG>& phi, // DG1 (3)
    const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie)
{
    if (smesh.landmask[c1] == 0)
        return;
    if (smesh.landmask[c2] == 0)
        return;

    const LocalEdgeVector<GAUSSPOINTS1D(DG)> vel_gauss
        = normalvel_Y.row(ie) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>;

    const LocalEdgeVector<GAUSSPOINTS1D(DG)> tmp = (vel_gauss.cwiseMax(0).array()
            * (rightedgeofcell<DG>(phi, c1) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>).array()
        + vel_gauss.cwiseMin(0).array()
            * (leftedgeofcell<DG>(phi, c2) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>).array());

    // - [[psi]] sind we're on the left side
    phiup.row(c1) -= dt * tmp * PSIe_w<DG, GAUSSPOINTS1D(DG), 1>;
    phiup.row(c2) += dt * tmp * PSIe_w<DG, GAUSSPOINTS1D(DG), 3>;
}

// Cell term with the full velocity: same matrices as
// ParametricTransportMap::InitializeAdvectionCellTerms, with NGP_CELL Gauss points
template <int DG> void DGTransport<DG>::InitializeCellTermsFull()
{
    constexpr int Q = NGP_CELL;
    cellTermFullX.resize(smesh.nelements);
    cellTermFullY.resize(smesh.nelements);

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        const Eigen::Matrix<FloatType, 2, Q * Q> dxT
            = ParametricTools::dxT<Q>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<Q>.array();
        const Eigen::Matrix<FloatType, 2, Q * Q> dyT
            = ParametricTools::dyT<Q>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<Q>.array();

        // [J dT^{-T} nabla phi]_1, [J dT^{-T} nabla phi]_2
        cellTermFullX[eid] = PSIx<DG, Q>.array().rowwise() * dyT.row(1).array()
            - PSIy<DG, Q>.array().rowwise() * dxT.row(1).array();
        cellTermFullY[eid] = PSIy<DG, Q>.array().rowwise() * dxT.row(0).array()
            - PSIx<DG, Q>.array().rowwise() * dyT.row(0).array();
    }
}

template <int DG>
void DGTransport<DG>::cell_term_full(
    DGVector<DG>& phiup, const DGVector<DG>& phi, const FloatType dt, const size_t eid)
{
    if constexpr (DG == 1) {
        return; // dG0: grad psi = 0
    } else {
        if (smesh.landmask[eid] == 0)
            return;
        constexpr int Q = NGP_CELL;
        const Eigen::Matrix<FloatType, 1, Q * Q> phi_gauss = phi.row(eid) * PSI<DG, Q>;
        phiup.row(eid) += dt
            * (cellTermFullX[eid].array().rowwise() * cellVel_X.row(eid).array()
                + cellTermFullY[eid].array().rowwise() * cellVel_Y.row(eid).array())
                  .matrix()
            * phi_gauss.transpose();
    }
}

// Edge terms with the full velocity: u_hat.n * length is given in the NGP_EDGE Gauss
// points, upwind in each Gauss point
template <int DG>
void DGTransport<DG>::edge_term_X_full(DGVector<DG>& phiup, const DGVector<DG>& phi,
    const FloatType dt, const size_t c1, const size_t c2, const size_t ie)
{
    if (smesh.landmask[c1] == 0 || smesh.landmask[c2] == 0)
        return;
    constexpr int NGP = NGP_EDGE;
    const LocalEdgeVector<NGP> vel = edgeNormalVel_X.row(ie);
    // phi of c1 (below) where u_hat.n > 0, of c2 (above) where < 0
    const LocalEdgeVector<NGP> flux = (vel.array().max(0)
            * (topedgeofcell<DG>(phi, c1) * PSIe<EDGEDOFS(DG), NGP>).array()
        + vel.array().min(0) * (bottomedgeofcell<DG>(phi, c2) * PSIe<EDGEDOFS(DG), NGP>).array())
                                          .matrix();
    phiup.row(c1) -= dt * flux * PSIe_w<DG, NGP, 2>;
    phiup.row(c2) += dt * flux * PSIe_w<DG, NGP, 0>;
}

template <int DG>
void DGTransport<DG>::edge_term_Y_full(DGVector<DG>& phiup, const DGVector<DG>& phi,
    const FloatType dt, const size_t c1, const size_t c2, const size_t ie)
{
    if (smesh.landmask[c1] == 0 || smesh.landmask[c2] == 0)
        return;
    constexpr int NGP = NGP_EDGE;
    const LocalEdgeVector<NGP> vel = edgeNormalVel_Y.row(ie);
    // phi of c1 (left) where u_hat.n > 0, of c2 (right) where < 0
    const LocalEdgeVector<NGP> flux
        = (vel.array().max(0) * (rightedgeofcell<DG>(phi, c1) * PSIe<EDGEDOFS(DG), NGP>).array()
            + vel.array().min(0) * (leftedgeofcell<DG>(phi, c2) * PSIe<EDGEDOFS(DG), NGP>).array())
              .matrix();
    phiup.row(c1) -= dt * flux * PSIe_w<DG, NGP, 1>;
    phiup.row(c2) += dt * flux * PSIe_w<DG, NGP, 3>;
}

// Outer boundary edge: upwind flux as in edge_term_X/Y_full, the missing neighbour is the
// exterior state boundaryInflowValue. Normals to +y / +x as in the interior.
template <int DG>
template <int Side>
void DGTransport<DG>::boundary_term_full(
    DGVector<DG>& phiup, const DGVector<DG>& phi, const FloatType dt, const size_t cell)
{
    if (smesh.landmask[cell] == 0)
        return;
    constexpr int NGP = NGP_EDGE;
    const size_t nx = smesh.nx;
    const size_t ix = cell % nx;
    const size_t iy = cell / nx;
    const LocalEdgeVector<NGP> phiExt = LocalEdgeVector<NGP>::Constant(boundaryInflowValue);

    if constexpr (Side == 0) { // bottom: exterior below, cell above
        const LocalEdgeVector<NGP> vel = edgeNormalVel_X.row(cell);
        const LocalEdgeVector<NGP> phiCell
            = bottomedgeofcell<DG>(phi, cell) * PSIe<EDGEDOFS(DG), NGP>;
        const LocalEdgeVector<NGP> flux
            = (vel.array().max(0) * phiExt.array() + vel.array().min(0) * phiCell.array()).matrix();
        phiup.row(cell) += dt * flux * PSIe_w<DG, NGP, 0>;
    } else if constexpr (Side == 1) { // right: cell left, exterior right
        const LocalEdgeVector<NGP> vel = edgeNormalVel_Y.row(iy * (nx + 1) + ix + 1);
        const LocalEdgeVector<NGP> phiCell
            = rightedgeofcell<DG>(phi, cell) * PSIe<EDGEDOFS(DG), NGP>;
        const LocalEdgeVector<NGP> flux
            = (vel.array().max(0) * phiCell.array() + vel.array().min(0) * phiExt.array()).matrix();
        phiup.row(cell) -= dt * flux * PSIe_w<DG, NGP, 1>;
    } else if constexpr (Side == 2) { // top: cell below, exterior above
        const LocalEdgeVector<NGP> vel = edgeNormalVel_X.row(cell + nx);
        const LocalEdgeVector<NGP> phiCell = topedgeofcell<DG>(phi, cell) * PSIe<EDGEDOFS(DG), NGP>;
        const LocalEdgeVector<NGP> flux
            = (vel.array().max(0) * phiCell.array() + vel.array().min(0) * phiExt.array()).matrix();
        phiup.row(cell) -= dt * flux * PSIe_w<DG, NGP, 2>;
    } else { // left: exterior left, cell right
        const LocalEdgeVector<NGP> vel = edgeNormalVel_Y.row(iy * (nx + 1) + ix);
        const LocalEdgeVector<NGP> phiCell
            = leftedgeofcell<DG>(phi, cell) * PSIe<EDGEDOFS(DG), NGP>;
        const LocalEdgeVector<NGP> flux
            = (vel.array().max(0) * phiExt.array() + vel.array().min(0) * phiCell.array()).matrix();
        phiup.row(cell) += dt * flux * PSIe_w<DG, NGP, 3>;
    }
}

template <int DG>
void DGTransport<DG>::DGTransportOperator(const ParametricMesh& smesh, const FloatType dt,
    const DGVector<DG>& vx, const DGVector<DG>& vy, const EdgeVector<EDGEDOFS(DG)>& normalvel_X,
    const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const DGVector<DG>& phi, DGVector<DG>& phiup)
{
    phiup.zero();

    // Cell terms
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        if (projectVelocity)
            cell_term(smesh, dt, phiup, phi, vx, vy, eid);
        else
            cell_term_full(phiup, phi, dt, eid);

    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        size_t ic = iy * smesh.nx; // first index of left cell in row
        size_t ie = iy * (smesh.nx + 1) + 1; // first index of inner velocity in row
        for (size_t i = 0; i < smesh.nx - 1; ++i, ++ic, ++ie)
            if (projectVelocity)
                edge_term_Y(smesh, dt, phiup, phi, normalvel_Y, ic, ic + 1, ie);
            else
                edge_term_Y_full(phiup, phi, dt, ic, ic + 1, ie);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < smesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        size_t ie = ix + smesh.nx; // first index of inner velocity in column
        for (size_t i = 0; i < smesh.ny - 1; ++i, ic += smesh.nx, ie += smesh.nx)
            if (projectVelocity)
                edge_term_X(smesh, dt, phiup, phi, normalvel_X, ic, ic + smesh.nx, ie);
            else
                edge_term_X_full(phiup, phi, dt, ic, ic + smesh.nx, ie);
    }

    // Periodic
    for (size_t pc = 0; pc < smesh.periodic.size(); ++pc) {
#pragma omp parallel for
        for (size_t i = 0; i < smesh.periodic[pc].size(); ++i) {
            if (smesh.periodic[pc][i][0] == 0) // X-edge (bottom / top)
                edge_term_X(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc][i][1],
                    smesh.periodic[pc][i][2], smesh.periodic[pc][i][3]);
            else if (smesh.periodic[pc][i][0] == 1) // Y-edge (left / right)
                edge_term_Y(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc][i][1],
                    smesh.periodic[pc][i][2], smesh.periodic[pc][i][3]);
            else {
                std::cerr << "Wrong periodic boundary information in the mesh. Boundary side "
                          << smesh.periodic[pc][i][0] << " not valid" << std::endl;
                abort();
            }
        }
    }

    // TR 07.04.2025
    //
    // Dirichlet: we do not need to do anything special in DGTransport. On Dirichlet
    // boundaries, the normal velocity is zero. Hence, any flux is zero.

    //
    // Neumann inflow / outflow
    // This can only appear on outer mesh-boundaries, whenever the element at the edge
    // is ice.

    // TO DISCUSS: there can't be in/outflow and periodic on the same edge. We either need a
    // structure for in/outflow or we say: either there is an in/outflow boundary or periodic. But
    // not both..
    if (!projectVelocity) {
        // full velocity: upwind with the exterior state on the outer boundary
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nx; ++i) {
            boundary_term_full<0>(phiup, phi, dt, i);
            boundary_term_full<2>(phiup, phi, dt, i + smesh.nx * (smesh.ny - 1));
        }
#pragma omp parallel for
        for (size_t i = 0; i < smesh.ny; ++i) {
            boundary_term_full<3>(phiup, phi, dt, i * smesh.nx);
            boundary_term_full<1>(phiup, phi, dt, (i + 1) * smesh.nx - 1);
        }
    } else if (smesh.periodic.size() == 0) {
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nx; ++i) {
            const size_t element_lower = i;
            const size_t edge_lower = i;
            const size_t element_upper = i + smesh.nx * (smesh.ny - 1);
            const size_t edge_upper = i + smesh.nx * smesh.ny;
            if (smesh.landmask[element_lower] == 1)
                boundary_lower(smesh, dt, phiup, phi, normalvel_X, element_lower, edge_lower);
            if (smesh.landmask[element_upper] == 1)
                boundary_upper(smesh, dt, phiup, phi, normalvel_X, element_upper, edge_upper);
        }
#pragma omp parallel for
        for (size_t i = 0; i < smesh.ny; ++i) {
            const size_t element_left = i * smesh.nx;
            const size_t edge_left = i * (smesh.nx + 1);
            const size_t element_right = (i + 1) * smesh.nx - 1;
            const size_t edge_right = i * (smesh.nx + 1) + smesh.nx;
            if (smesh.landmask[element_left] == 1)
                boundary_left(smesh, dt, phiup, phi, normalvel_Y, element_left, edge_left);
            if (smesh.landmask[element_right] == 1)
                boundary_right(smesh, dt, phiup, phi, normalvel_Y, element_right, edge_right);
        }
    }

#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        phiup.row(eid) = parammap.InverseDGMassMatrix[eid] * phiup.row(eid).transpose();
}

template <int DG> void DGTransport<DG>::step_rk1(const FloatType dt, DGVector<DG>& phi)
{
    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1);

    phi += tmp1;
}

template <int DG> void DGTransport<DG>::step_rk2(const FloatType dt, DGVector<DG>& phi)
{
    DGTransportOperator(
        smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1); // tmp1 = k * F(u)

    phi += tmp1; // phi = phi + k * F(u)     (i.e.: implicit Euler)

    DGTransportOperator(
        smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp2); // tmp1 = k * F( u + k * F(u) )

    phi += 0.5_ft * (tmp2 - tmp1);
}

template <int DG> void DGTransport<DG>::step_rk3(const FloatType dt, DGVector<DG>& phi)
{
    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
        tmp1); // tmp1 = k * F(u)  // K1 in Heun(3)
    tmp1 += phi; // phi + h f(phi)

    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, tmp1,
        tmp2); // tmp2 = f( u + h f(u) )
    tmp2 += tmp1;
    tmp2 *= 0.25;
    tmp2 += 0.75 * phi;

    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, tmp2,
        tmp3); // k * F(k1) // K2 in Heun(3)
    tmp3 += tmp2;

    phi *= 1.0 / 3.0;
    phi += 2.0 / 3.0 * tmp3;
}

template <int DG> void DGTransport<DG>::step(const FloatType dt, DGVector<DG>& phi)
{
    if (timesteppingscheme == "rk1")
        step_rk1(dt, phi);
    else if (timesteppingscheme == "rk2")
        step_rk2(dt, phi);
    else if (timesteppingscheme == "rk3")
        step_rk3(dt, phi);
    else {
        std::cerr << "Time stepping scheme '" << timesteppingscheme << "' not known!" << std::endl;
        abort();
    }
}

template class DGTransport<1>;
template class DGTransport<3>;
template class DGTransport<4>;
template class DGTransport<6>;
template class DGTransport<8>;

template void DGTransport<1>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<1>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<3>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<3>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<4>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<4>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<6>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<6>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<8>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);

// transport with the DG velocity of the LDG momentum
template void DGTransport<1>::prepareAdvectionDG(
    const DGVector<3>& dg_vx, const DGVector<3>& dg_vy, const FloatType a);
template void DGTransport<1>::prepareAdvectionDG(
    const DGVector<4>& dg_vx, const DGVector<4>& dg_vy, const FloatType a);
template void DGTransport<1>::prepareAdvectionDG(
    const DGVector<6>& dg_vx, const DGVector<6>& dg_vy, const FloatType a);
template void DGTransport<1>::prepareAdvectionDG(
    const DGVector<9>& dg_vx, const DGVector<9>& dg_vy, const FloatType a);
template void DGTransport<3>::prepareAdvectionDG(
    const DGVector<3>& dg_vx, const DGVector<3>& dg_vy, const FloatType a);
template void DGTransport<3>::prepareAdvectionDG(
    const DGVector<4>& dg_vx, const DGVector<4>& dg_vy, const FloatType a);
template void DGTransport<3>::prepareAdvectionDG(
    const DGVector<6>& dg_vx, const DGVector<6>& dg_vy, const FloatType a);
template void DGTransport<3>::prepareAdvectionDG(
    const DGVector<9>& dg_vx, const DGVector<9>& dg_vy, const FloatType a);
template void DGTransport<4>::prepareAdvectionDG(
    const DGVector<3>& dg_vx, const DGVector<3>& dg_vy, const FloatType a);
template void DGTransport<4>::prepareAdvectionDG(
    const DGVector<4>& dg_vx, const DGVector<4>& dg_vy, const FloatType a);
template void DGTransport<4>::prepareAdvectionDG(
    const DGVector<6>& dg_vx, const DGVector<6>& dg_vy, const FloatType a);
template void DGTransport<4>::prepareAdvectionDG(
    const DGVector<9>& dg_vx, const DGVector<9>& dg_vy, const FloatType a);
template void DGTransport<6>::prepareAdvectionDG(
    const DGVector<3>& dg_vx, const DGVector<3>& dg_vy, const FloatType a);
template void DGTransport<6>::prepareAdvectionDG(
    const DGVector<4>& dg_vx, const DGVector<4>& dg_vy, const FloatType a);
template void DGTransport<6>::prepareAdvectionDG(
    const DGVector<6>& dg_vx, const DGVector<6>& dg_vy, const FloatType a);
template void DGTransport<6>::prepareAdvectionDG(
    const DGVector<9>& dg_vx, const DGVector<9>& dg_vy, const FloatType a);

} /* namespace Nextsim */
