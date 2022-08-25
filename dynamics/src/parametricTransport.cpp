/*!
 * @file ParametricTransport.cpp
 * @date July 10, 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "ParametricTransport.hpp"
#include "ParametricTools.hpp"
//#include "dgTimeStepping.hpp"
#include "stopwatch.hpp"

namespace Nextsim {

extern Timer GlobalTimer;

//! returns the localization of the cell vector to the edges
template <int DGcell, int DGedge>
Eigen::Matrix<Nextsim::FloatType, 1, DGedge>
leftedgeofcell(const CellVector<DGcell>& cv, size_t eid);
template <int DGcell, int DGedge>
Eigen::Matrix<Nextsim::FloatType, 1, DGedge>
rightedgeofcell(const CellVector<DGcell>& cv, size_t eid);
template <int DGcell, int DGedge>
Eigen::Matrix<Nextsim::FloatType, 1, DGedge>
bottomedgeofcell(const CellVector<DGcell>& cv, size_t eid);
template <int DGcell, int DGedge>
Eigen::Matrix<Nextsim::FloatType, 1, DGedge>
topedgeofcell(const CellVector<DGcell>& cv, size_t eid);

// dG0 (1 in cell, 1 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
leftedgeofcell(const CellVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
rightedgeofcell(const CellVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
bottomedgeofcell(const CellVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
topedgeofcell(const CellVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}

// dG1 (3 in cell, 2 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
leftedgeofcell(const CellVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
rightedgeofcell(const CellVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
bottomedgeofcell(const CellVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
topedgeofcell(const CellVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
leftedgeofcell(const CellVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5), cv(eid, 4));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
rightedgeofcell(const CellVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5), cv(eid, 4));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
bottomedgeofcell(const CellVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5), cv(eid, 3));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
topedgeofcell(const CellVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5), cv(eid, 3));
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <int DGcell, int DGedge>
void ParametricTransport<DGcell, DGedge>::reinitnormalvelocity()
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
            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_left = smesh.edgevector(ey, ey + smesh.nx + 1);

            normalvel_Y.row(ey) += 0.5 * (tangent_left(0, 1) * leftedgeofcell<DGcell, DGedge>(velx, cy) - tangent_left(0, 0) * leftedgeofcell<DGcell, DGedge>(vely, cy));

            // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_right = smesh.edgevector(ey + 1, ey + smesh.nx + 2);

            normalvel_Y.row(ey + 1) += 0.5 * (tangent_right(0, 1) * rightedgeofcell<DGcell, DGedge>(velx, cy) - tangent_right(0, 0) * rightedgeofcell<DGcell, DGedge>(vely, cy));
        }

        // scale boundary
        normalvel_Y.row(iy * (smesh.nx + 1)) *= 2.0;
        normalvel_Y.row((iy + 1) * (smesh.nx + 1) - 1) *= 2.0;
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
            // un-normed tangent vector of bottom edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_bottom = smesh.edgevector(nx, nx + 1);

            normalvel_X.row(cx) += 0.5 * (-tangent_bottom(0, 1) * bottomedgeofcell<DGcell, DGedge>(velx, cx) + tangent_bottom(0, 0) * bottomedgeofcell<DGcell, DGedge>(vely, cx));

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_top = smesh.edgevector(nx + smesh.nx + 1, nx + smesh.nx + 2);

            normalvel_X.row(cx + smesh.nx) += 0.5 * (-tangent_top(0, 1) * topedgeofcell<DGcell, DGedge>(velx, cx) + tangent_top(0, 0) * topedgeofcell<DGcell, DGedge>(vely, cx));
        }

        // scale boundary
        normalvel_X.row(ix) *= 2.0;
        normalvel_X.row(ix + smesh.ny * smesh.nx) *= 2.0;
    }
}

////////////////////////////////////////////////// CELL TERM

template <int DG>
void cell_term(const SasipMesh& smesh, double dt,
    CellVector<DG>& phiup, const CellVector<DG>& phi,
    const CellVector<DG>& vx,
    const CellVector<DG>& vy, const size_t ic);
template <>
void cell_term(const SasipMesh& smesh, double dt,
    CellVector<1>& phiup, const CellVector<1>& phi,
    const CellVector<1>& vx,
    const CellVector<1>& vy, const size_t ic) { }

// higher order terms with gauss quadrature
//
//     // - (Psi v, \nabla phi)
// -  wq *  Psi[q] * v[q] * ( Jq * JT^{-T} ) [BIGx, BIGy]

template <int DG>
void cell_term(const SasipMesh& smesh, double dt,
    CellVector<DG>& phiup,
    const CellVector<DG>& phi,
    const CellVector<DG>& vx,
    const CellVector<DG>& vy, const size_t eid)
{
#define NGP 3
    // - w_q Psi(q) * ( vx(q) dx_phi_i(q) + vy(q) dy_phi_i(q) )

    // REICHT vielleicht 2-Punkt Gauss?

    const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> vx_gauss = vx.row(eid) * PSI<DG, NGP>; //!< velocity in GP
    const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> vy_gauss = vy.row(eid) * PSI<DG, NGP>;

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
    const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> phi_gauss = GAUSSWEIGHTS<NGP>.array() * (phi.row(eid) * PSI<DG, NGP>).array();

    const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dxT = ParametricTools::dxT<NGP>(smesh, eid);
    const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dyT = ParametricTools::dyT<NGP>(smesh, eid);

    // [J dT^{-T} nabla phi]_1
    phiup.row(eid) += dt * ((PSIx<DG, NGP>.array().rowwise() * dyT.row(1).array() - PSIy<DG, NGP>.array().rowwise() * dxT.row(1).array()).rowwise() * vx_gauss.array()
                          // [J dT^{-T} nabla phi]_2
                          + (PSIy<DG, NGP>.array().rowwise() * dxT.row(0).array() - PSIx<DG, NGP>.array().rowwise() * dyT.row(0).array()).rowwise() * vy_gauss.array())
                               .matrix()
        * phi_gauss.transpose();

    //  - x * 1./8000

    //   // vx, vy dG(0) !!!
    //   // \nabla phi = (0, 1, 1)

    //   // AV[0] * dx Phi(1) * VX[0]
    //   phiup(ic, 1) += inversemasscell(1) / mesh.hx
    //       * (phi(ic, 0) * vx(ic, 0) + 1. / 12. * phi(ic, 1) * vx(ic, 1)
    //           + 1. / 12. * phi(ic, 2) * vx(ic, 2));
    //   // AV[0] * dy Phi(2) * VX[0]
    //   phiup(ic, 2) += inversemasscell(2) / mesh.hy
    //       * (phi(ic, 0) * vy(ic, 0) + 1. / 12. * phi(ic, 1) * vy(ic, 1)
    //           + 1. / 12. * phi(ic, 2) * vy(ic, 2));

    // alle anderen Beitraege 0

#undef NGP
}
////////////////////////////////////////////////// BOUNDARY HANDLING

void boundary_lower(const SasipMesh& smesh, const double dt, CellVector<1>& phiup,
    const CellVector<1>& phi, const EdgeVector<1>& normalvel_X, const size_t c,
    const size_t e)
{
    phiup(c, 0) -= dt * std::max(0., -normalvel_X(e, 0)) * phi(c, 0);
}
template <int DG, int DGedge>
void boundary_lower(const SasipMesh& smesh, const double dt, CellVector<DG>& phiup,
    const CellVector<DG>& phi, const EdgeVector<DGedge>& normalvel_X, const size_t c, const size_t e)
{
    // GP = DGEDGE
    LocalEdgeVector<DGedge> vel_gauss = normalvel_X.row(e) * PSIe<DGedge, DGedge>;
    // block<1, 2>(e, 0) * PSIe<2,2>;
    LocalEdgeVector<DGedge> tmp = ((bottomedgeofcell<DG, DGedge>(phi, c) * PSIe<DGedge, DGedge>).array() * (-vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, DGedge, 0>;
}

void boundary_upper(const SasipMesh& smesh, const double dt, CellVector<1>& phiup,
    const CellVector<1>& phi, const EdgeVector<1>& normalvel_X, const size_t c,
    const size_t e)
{
    phiup(c, 0) -= dt * std::max(0., normalvel_X(e, 0)) * phi(c, 0);
}
template <int DG, int DGedge>
void boundary_upper(const SasipMesh& smesh, const double dt, CellVector<DG>& phiup,
    const CellVector<DG>& phi, const EdgeVector<DGedge>& normalvel_X, const size_t c, const size_t e)
{
    LocalEdgeVector<DGedge> vel_gauss = normalvel_X.row(e) * PSIe<DGedge, DGedge>;
    LocalEdgeVector<DGedge> tmp = ((topedgeofcell<DG, DGedge>(phi, c) * PSIe<DGedge, DGedge>).array() * (vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, DGedge, 2>;
}

void boundary_left(const SasipMesh& smesh, const double dt, CellVector<1>& phiup,
    const CellVector<1>& phi, const EdgeVector<1>& normalvel_Y, const size_t c,
    const size_t e)
{
    phiup(c, 0) -= dt * std::max(0., -normalvel_Y(e, 0)) * phi(c, 0);
}
template <int DG, int DGedge>
void boundary_left(const SasipMesh& smesh, const double dt, CellVector<DG>& phiup,
    const CellVector<DG>& phi, const EdgeVector<DGedge>& normalvel_Y, const size_t c, const size_t e)
{
    LocalEdgeVector<DGedge> vel_gauss = normalvel_Y.row(e) * PSIe<DGedge, DGedge>;
    LocalEdgeVector<DGedge> tmp = ((leftedgeofcell<DG, DGedge>(phi, c) * PSIe<DGedge, DGedge>).array() * (-vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, DGedge, 3>;
}
void boundary_right(const SasipMesh& smesh, const double dt, CellVector<1>& phiup,
    const CellVector<1>& phi, const EdgeVector<1>& normalvel_Y, const size_t c,
    const size_t e)
{
    phiup(c, 0) -= dt * std::max(0., normalvel_Y(e, 0)) * phi(c, 0);
}
template <int DG, int DGedge>
void boundary_right(const SasipMesh& smesh, const double dt, CellVector<DG>& phiup,
    const CellVector<DG>& phi, const EdgeVector<DGedge>& normalvel_Y, const size_t c, const size_t e)
{
    LocalEdgeVector<DGedge> vel_gauss = normalvel_Y.row(e) * PSIe<DGedge, DGedge>;
    LocalEdgeVector<DGedge> tmp = ((rightedgeofcell<DG, DGedge>(phi, c) * PSIe<DGedge, DGedge>).array() * (vel_gauss.array().max(0)));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, DGedge, 1>;
}

void edge_term_X(const SasipMesh& smesh, const double dt, CellVector<1>& phiup, const CellVector<1>& phi, // DG0 (1)
    const EdgeVector<1>& normalvel_X, const size_t c1, const size_t c2, const size_t ie)
{
    double bottom = phi(c1, 0);
    double top = phi(c2, 0);
    double vel = normalvel_X(ie, 0);

    phiup(c1, 0) -= dt * (std::max(vel, 0.) * bottom + std::min(vel, 0.) * top);
    phiup(c2, 0) += dt * (std::max(vel, 0.) * bottom + std::min(vel, 0.) * top);
}
void edge_term_Y(const SasipMesh& smesh, const double dt, CellVector<1>& phiup, const CellVector<1>& phi, // DG0 (1)
    const EdgeVector<1>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie)
{
    double left = phi(c1, 0);
    double right = phi(c2, 0);
    double vel = normalvel_Y(ie, 0);

    phiup(c1, 0) -= dt * (std::max(vel, 0.) * left + std::min(vel, 0.) * right);
    phiup(c2, 0) += dt * (std::max(vel, 0.) * left + std::min(vel, 0.) * right);
}

template <int DG, int DGedge>
void edge_term_X(const SasipMesh& smesh, const double dt, CellVector<DG>& phiup, const CellVector<DG>& phi, // DG1 (3)
    const EdgeVector<DGedge>& normalvel_X, const size_t c1, const size_t c2, const size_t ie)
{
    const LocalEdgeVector<DGedge> vel_gauss = normalvel_X.row(ie) * PSIe<DGedge, DGedge>;

    const LocalEdgeVector<DGedge> tmp = (vel_gauss.array().max(0) * (topedgeofcell<DG, DGedge>(phi, c1) * PSIe<DGedge, DGedge>).array()
        + vel_gauss.array().min(0) * (bottomedgeofcell<DG, DGedge>(phi, c2) * PSIe<DGedge, DGedge>).array());
    phiup.row(c1) -= dt * tmp * PSIe_w<DG, DGedge, 2>;
    phiup.row(c2) += dt * tmp * PSIe_w<DG, DGedge, 0>;
}
template <int DG, int DGedge>
void edge_term_Y(const SasipMesh& smesh, const double dt, CellVector<DG>& phiup, const CellVector<DG>& phi, // DG1 (3)
    const EdgeVector<DGedge>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie)
{
    const LocalEdgeVector<DGedge> vel_gauss = normalvel_Y.row(ie) * PSIe<DGedge, DGedge>;
    const LocalEdgeVector<DGedge> tmp = (vel_gauss.array().max(0) * (rightedgeofcell<DG, DGedge>(phi, c1) * PSIe<DGedge, DGedge>).array()
        + vel_gauss.array().min(0) * (leftedgeofcell<DG, DGedge>(phi, c2) * PSIe<DGedge, DGedge>).array());

    // - [[psi]] sind we're on the left side
    phiup.row(c1) -= dt * tmp * PSIe_w<DG, DGedge, 1>;
    phiup.row(c2) += dt * tmp * PSIe_w<DG, DGedge, 3>;
}

template <int DGcell, int DGedge>
void parametricTransportOperator(const SasipMesh& smesh, const double dt,
    const CellVector<DGcell>& vx,
    const CellVector<DGcell>& vy,
    const EdgeVector<DGedge>& normalvel_X,
    const EdgeVector<DGedge>& normalvel_Y,
    const CellVector<DGcell>& phi, CellVector<DGcell>& phiup)
{
    phiup.zero();

    GlobalTimer.start("-- -- --> cell term");
    // Cell terms
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        cell_term<DGcell>(smesh, dt, phiup, phi, vx, vy, eid);
    GlobalTimer.stop("-- -- --> cell term");

    GlobalTimer.start("-- -- --> edge terms");
    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        size_t ic = iy * smesh.nx; // first index of left cell in row
        size_t ie = iy * (smesh.nx + 1) + 1; // first index of inner velocity in row

        for (size_t i = 0; i < smesh.nx - 1; ++i, ++ic, ++ie)
            edge_term_Y(smesh, dt, phiup, phi, normalvel_Y, ic, ic + 1, ie);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < smesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        size_t ie = ix + smesh.nx; // first index of inner velocity in column
        for (size_t i = 0; i < smesh.ny - 1; ++i, ic += smesh.nx, ie += smesh.nx)
            edge_term_X(smesh, dt, phiup, phi, normalvel_X, ic, ic + smesh.nx, ie);
    }
    GlobalTimer.stop("-- -- --> edge terms");

    // boundaries
    GlobalTimer.start("-- -- --> boundaries");
    // lower & upper
    //#pragma omp parallel for
    const size_t eupper0 = smesh.nx * smesh.ny;
    for (size_t ix = 0; ix < smesh.nx; ++ix) {
        const size_t clower = ix;
        const size_t elower = ix;
        const size_t cupper = smesh.nelements - smesh.nx + ix;
        const size_t eupper = eupper0 + ix;

        boundary_lower(smesh, dt, phiup, phi, normalvel_X, clower, elower);
        boundary_upper(smesh, dt, phiup, phi, normalvel_X, cupper, eupper);
    }
    // left & right
    //#pragma omp parallel for
    size_t eright0 = smesh.nx;

    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        const size_t cleft = iy * smesh.nx;
        const size_t eleft = iy * (smesh.nx + 1);
        const size_t cright = eright0 - 1 + iy * smesh.nx;
        const size_t eright = eright0 + iy * (smesh.nx + 1);

        boundary_left(smesh, dt, phiup, phi, normalvel_Y, cleft, eleft);
        boundary_right(smesh, dt, phiup, phi, normalvel_Y, cright, eright);
    }
    GlobalTimer.stop("-- -- --> boundaries");

    GlobalTimer.start("-- -- --> inverse mass");
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        phiup.row(eid) = (ParametricTools::massMatrix<DGcell>(smesh, eid).inverse() * phiup.row(eid).transpose());
    }

    GlobalTimer.stop("-- -- --> inverse mass");
}

template <int DGcell, int DGedge>
void ParametricTransport<DGcell, DGedge>::step_rk1(const double dt, CellVector<DGcell>& phi)
{
    parametricTransportOperator<DGcell, DGedge>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1);

    phi += tmp1;
}

template <int DGcell, int DGedge>
void ParametricTransport<DGcell, DGedge>::step_rk2(const double dt, CellVector<DGcell>& phi)
{
    parametricTransportOperator<DGcell, DGedge>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1); // tmp1 = k * F(u)

    phi += tmp1; // phi = phi + k * F(u)     (i.e.: implicit Euler)

    parametricTransportOperator<DGcell, DGedge>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp2); // tmp1 = k * F( u + k * F(u) )

    phi += 0.5 * (tmp2 - tmp1);
}

template <int DGcell, int DGedge>
void ParametricTransport<DGcell, DGedge>::step_rk3(const double dt, CellVector<DGcell>& phi)
{
    parametricTransportOperator<DGcell, DGedge>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
        tmp1); // tmp1 = k * F(u)  // K1 in Heun(3)

    phi += 1. / 3. * tmp1; // phi = phi + k/3 * F(u)   (i.e.: implicit Euler)
    parametricTransportOperator<DGcell, DGedge>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
        tmp2); // k * F(k1) // K2 in Heun(3)
    phi -= 1. / 3. * tmp1; // phi = phi + k/3 * F(u)   (i.e.: implicit Euler)

    phi += 2. / 3. * tmp2;
    parametricTransportOperator<DGcell, DGedge>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
        tmp3); // k * F(k2) // K3 in Heun(3)
    phi -= 2. / 3. * tmp2;

    phi += 0.25 * tmp1 + 0.75 * tmp3;
}

template <int DGcell, int DGedge>
void ParametricTransport<DGcell, DGedge>::step(const double dt, CellVector<DGcell>& phi)
{
    GlobalTimer.start("-- --> step");
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
    GlobalTimer.stop("-- --> step");
}

template class ParametricTransport<1, 1>;
template class ParametricTransport<3, 2>;
template class ParametricTransport<6, 3>;

} /* namespace Nextsim */
