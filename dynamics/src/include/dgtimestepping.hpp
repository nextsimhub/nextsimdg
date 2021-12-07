/*----------------------------   dgtimestepping.h ---------------------------*/
#ifndef __dgtimestepping_H
#define __dgtimestepping_H
/*----------------------------   dgtimestepping.h ---------------------------*/

#include "dgvector_manipulations.hpp"
#include "stopwatch.hpp"
#include "timemesh.hpp"

namespace Nextsim {

extern Timer GlobalTimer;

// pre-computed matrices for assembling dG-transport
// the Gauss rule for dG(n) is set to n+1 points
// this might not be enough?
#include "dgbasisfunctions_gausspoints.hpp"

// compute the cell terms starting dG(1)
template <int DGdegree>
void cell_term(const Mesh& mesh,
    const LocalCellVector<DGdegree> inversemasscell,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const CellVector<DGdegree>& vx,
    const CellVector<DGdegree>& vy,
    const size_t ic);
template <>
void cell_term(const Mesh& mesh,
    const LocalCellVector<0> inversemasscell,
    CellVector<0>& phiup,
    const CellVector<0>& phi,
    const CellVector<0>& vx,
    const CellVector<0>& vy,
    const size_t ic)
{
}
template <>
void cell_term(const Mesh& mesh,
    const LocalCellVector<1> inversemasscell,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const CellVector<1>& vx,
    const CellVector<1>& vy,
    const size_t ic)
{
    // - (Psi v, \nabla phi)

    // vx, vy dG(0) !!!
    // \nabla phi = (0, 1, 1)

    // AV[0] * dx Phi(1) * VX[0]
    phiup(ic, 1) += inversemasscell(1) / mesh.hx * (phi(ic, 0) * vx(ic, 0) + 1. / 12. * phi(ic, 1) * vx(ic, 1) + 1. / 12. * phi(ic, 2) * vx(ic, 2));
    // AV[0] * dy Phi(2) * VX[0]
    phiup(ic, 2) += inversemasscell(2) / mesh.hy * (phi(ic, 0) * vy(ic, 0) + 1. / 12. * phi(ic, 1) * vy(ic, 1) + 1. / 12. * phi(ic, 2) * vy(ic, 2));

    // alle anderen Beitraege 0
}
template <>
void cell_term(const Mesh& mesh,
    const LocalCellVector<2> inversemasscell,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const CellVector<2>& vx,
    const CellVector<2>& vy,
    const size_t ic)
{
    // - (Psi v, \nabla phi)

    // =

    // - (psi * vx * d_x phi)
    // - (psi * vy * d_y phi)

    phiup(ic, 1) += inversemasscell(1) / mesh.hx * (phi(ic, 0) * vx(ic, 0) + 1. / 12. * phi(ic, 1) * vx(ic, 1) + 1. / 12. * phi(ic, 2) * vx(ic, 2) + 1. / 180. * phi(ic, 3) * vx(ic, 3) + 1. / 180. * phi(ic, 4) * vx(ic, 4) + 1. / 144. * phi(ic, 5) * vx(ic, 5));

    phiup(ic, 2) += inversemasscell(2) / mesh.hy * (phi(ic, 0) * vy(ic, 0) + 1. / 12. * phi(ic, 1) * vy(ic, 1) + 1. / 12. * phi(ic, 2) * vy(ic, 2) + 1. / 180. * phi(ic, 3) * vy(ic, 3) + 1. / 180. * phi(ic, 4) * vy(ic, 4) + 1. / 144. * phi(ic, 5) * vy(ic, 5));

    phiup(ic, 3) += inversemasscell(3) / mesh.hx * (1. / 6. * (phi(ic, 1) * vx(ic, 0) + phi(ic, 0) * vx(ic, 1)) + 1. / 90. * (phi(ic, 3) * vx(ic, 1) + phi(ic, 1) * vx(ic, 3)) + 1. / 72. * (phi(ic, 2) * vx(ic, 5) + phi(ic, 5) * vx(ic, 2)));

    phiup(ic, 4) += inversemasscell(4) / mesh.hy * (1. / 6. * (phi(ic, 2) * vy(ic, 0) + phi(ic, 0) * vy(ic, 2)) + 1. / 90. * (phi(ic, 2) * vy(ic, 4) + phi(ic, 4) * vy(ic, 2)) + 1. / 72. * (phi(ic, 1) * vy(ic, 5) + phi(ic, 5) * vy(ic, 1)));

    phiup(ic, 5) += inversemasscell(5) / mesh.hx * (1. / 12. * (phi(ic, 0) * vx(ic, 2) + phi(ic, 2) * vx(ic, 0)) + 1. / 144. * (phi(ic, 1) * vx(ic, 5) + phi(ic, 5) * vx(ic, 1)) + 1. / 180. * (phi(ic, 2) * vx(ic, 4) + phi(ic, 4) * vx(ic, 2)));

    phiup(ic, 5) += inversemasscell(5) / mesh.hy * (1. / 12. * (phi(ic, 0) * vy(ic, 1) + phi(ic, 1) * vy(ic, 0)) + 1. / 144. * (phi(ic, 2) * vy(ic, 5) + phi(ic, 5) * vy(ic, 2)) + 1. / 180. * (phi(ic, 1) * vy(ic, 3) + phi(ic, 3) * vy(ic, 1)));
    // all further contributions are zero
}

// compute the edge terms for the horizontal edges:  n = (0, +/- 1)
template <int DGdegree>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evy,
    const size_t c1,
    const size_t c2,
    const size_t e);

template <>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<0>& phiup,
    const CellVector<0>& phi,
    const EdgeVector<0>& evy,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    double bottom = phi(c1, 0);
    double top = phi(c2, 0);
    double vel = evy(ie, 0);

    phiup(c1, 0) -= timemesh.dt / mesh.hy * (std::max(vel, 0.) * bottom + std::min(vel, 0.) * top);
    phiup(c2, 0) += timemesh.dt / mesh.hy * (std::max(vel, 0.) * bottom + std::min(vel, 0.) * top);
}

template <>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evy,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    const LocalEdgeVector<1> bottom(
        phi(c1, 0) + 0.5 * phi(c1, 2),
        phi(c1, 1));
    const LocalEdgeVector<1> top(
        phi(c2, 0) - 0.5 * phi(c2, 2),
        phi(c2, 1));

    const LocalEdgeVector<1> vel_gauss = evy.block<1, 2>(ie, 0) * BiGe12;
    const LocalEdgeVector<1> tmp = (vel_gauss.array().max(0) * (bottom * BiGe12).array() + vel_gauss.array().min(0) * (top * BiGe12).array());
    phiup.block<1, 3>(c1, 0) -= timemesh.dt / mesh.hy * tmp * BiG12_2;
    phiup.block<1, 3>(c2, 0) += timemesh.dt / mesh.hy * tmp * BiG12_0;
}

template <>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evy,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    const LocalEdgeVector<2> bottom(
        phi(c1, 0) + 0.5 * phi(c1, 2) + 1. / 6. * phi(c1, 4),
        phi(c1, 1) + 0.5 * phi(c1, 5),
        phi(c1, 3));
    const LocalEdgeVector<2> top(
        phi(c2, 0) - 0.5 * phi(c2, 2) + 1. / 6. * phi(c2, 4),
        phi(c2, 1) - 0.5 * phi(c2, 5), phi(c2, 3));

    const LocalEdgeVector<2> vel_gauss = evy.block<1, 3>(ie, 0) * BiGe23;
    const LocalEdgeVector<2> tmp = (vel_gauss.array().max(0) * (bottom * BiGe23).array() + vel_gauss.array().min(0) * (top * BiGe23).array());
    phiup.block<1, 6>(c1, 0) -= timemesh.dt / mesh.hy * tmp * BiG23_2;
    phiup.block<1, 6>(c2, 0) += timemesh.dt / mesh.hy * tmp * BiG23_0;
}

// compute the edge terms for the vertical edges:  n = (+/- 1, 0)
template <int DGdegree>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<DGdegree>& phiup,

    const CellVector<DGdegree>& phi,

    const EdgeVector<DGdegree>& evx,
    const size_t c1,
    const size_t c2,
    const size_t e);

template <>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<0>& phiup,

    const CellVector<0>& phi,

    const EdgeVector<0>& evx,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    double left = phi(c1, 0);
    double right = phi(c2, 0);
    double vel = evx(ie, 0);

    phiup(c1, 0) -= timemesh.dt / mesh.hy * (std::max(vel, 0.) * left + std::min(vel, 0.) * right);
    phiup(c2, 0) += timemesh.dt / mesh.hy * (std::max(vel, 0.) * left + std::min(vel, 0.) * right);
}

template <>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evx,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    // average. cell: (1, x/h-1/2, y/h-1/2)  edge (1, x/h-1/2)
    const LocalEdgeVector<1> left(phi(c1, 0) + 0.5 * phi(c1, 1),
        phi(c1, 2));

    const LocalEdgeVector<1> right(phi(c2, 0) - 0.5 * phi(c2, 1),
        phi(c2, 2));

    const LocalEdgeVector<1> vel_gauss = evx.block<1, 2>(ie, 0) * BiGe12;
    const LocalEdgeVector<1> tmp = (vel_gauss.array().max(0) * (left * BiGe12).array() + vel_gauss.array().min(0) * (right * BiGe12).array());

    // - [[psi]] sind we're on the left side
    phiup.block<1, 3>(c1, 0) -= timemesh.dt / mesh.hx * tmp * BiG12_1;
    phiup.block<1, 3>(c2, 0) += timemesh.dt / mesh.hx * tmp * BiG12_3;
}

template <>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evx,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    // average. cell: (1, x/h-1/2, y/h-1/2)  edge (1, x/h-1/2)
    const LocalEdgeVector<2> left(phi(c1, 0) + 0.5 * phi(c1, 1) + 1. / 6. * phi(c1, 3),
        phi(c1, 2) + 0.5 * phi(c1, 5),
        phi(c1, 4));

    const LocalEdgeVector<2> right(phi(c2, 0) - 0.5 * phi(c2, 1) + 1. / 6. * phi(c2, 3),
        phi(c2, 2) - 0.5 * phi(c2, 5),
        phi(c2, 4));

    const LocalEdgeVector<2> vel_gauss = evx.block<1, 3>(ie, 0) * BiGe23;
    const LocalEdgeVector<2> tmp = (vel_gauss.array().max(0) * (left * BiGe23).array() + vel_gauss.array().min(0) * (right * BiGe23).array());

    // - [[psi]] sind we're on the left side
    phiup.block<1, 6>(c1, 0) -= timemesh.dt / mesh.hx * tmp * BiG23_1;
    phiup.block<1, 6>(c2, 0) += timemesh.dt / mesh.hx * tmp * BiG23_3;
}

//////////////////////////////////////////////////

//  <<  <v*n>^+ Psi Phi >>

template <int DGdegree>
void boundary_lower(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evy,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= timemesh.dt / mesh.hy * std::max(0., -evy(e, 0)) * phi(c, 0);
}
template <>
void boundary_lower(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evy,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<1> phi_lower
        = LocalEdgeVector<1>(phi(c, 0) - 0.5 * phi(c, 2), phi(c, 1)) * BiGe12;
    LocalEdgeVector<1> vel_gauss = evy.block<1, 2>(e, 0) * BiGe12;
    LocalEdgeVector<1> tmp = (phi_lower.array() * (-vel_gauss.array()).max(0));
    phiup.block<1, 3>(c, 0) -= timemesh.dt / mesh.hy * tmp * BiG12_0;
}
template <>
void boundary_lower(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evy,
    const size_t c,
    const size_t e)
{
    // 1, x-1/2, y-1/2, (x-1/2)^2-1/12, (y-1/2)^2-1/12, (x-1/2)(y-1/2)
    LocalEdgeVector<2> phi_lower
        = LocalEdgeVector<2>(phi(c, 0) - 0.5 * phi(c, 2) + 1. / 6. * phi(c, 4),
              phi(c, 1) - 0.5 * phi(c, 5),
              phi(c, 3))
        * BiGe23;
    // LocalEdgeVector<2> phi_gauss = phi_lower * BiGe23; // X
    LocalEdgeVector<2> vel_gauss = evy.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_lower.array() * (-vel_gauss.array()).max(0));
    phiup.block<1, 6>(c, 0) -= timemesh.dt / mesh.hy * tmp * BiG23_0;
}
template <int DGdegree>
void boundary_upper(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evy,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= timemesh.dt / mesh.hy * std::max(0., evy(e, 0)) * phi(c, 0);
}
template <>
void boundary_upper(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evy,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<1> phi_upper
        = LocalEdgeVector<1>(phi(c, 0) + 0.5 * phi(c, 2),
              phi(c, 1))
        * BiGe12;

    // LocalEdgeVector<2> phi_gauss = phi_upper * BiGe23; // X
    LocalEdgeVector<1> vel_gauss = evy.block<1, 2>(e, 0) * BiGe12;
    LocalEdgeVector<1> tmp = (phi_upper.array() * (vel_gauss.array()).max(0));
    phiup.block<1, 3>(c, 0) -= timemesh.dt / mesh.hy * tmp * BiG12_2;
}
template <>
void boundary_upper(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evy,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<2> phi_upper
        = LocalEdgeVector<2>(phi(c, 0) + 0.5 * phi(c, 2) + 1. / 6. * phi(c, 4),
              phi(c, 1) + 0.5 * phi(c, 5),
              phi(c, 3))
        * BiGe23;

    // LocalEdgeVector<2> phi_gauss = phi_upper * BiGe23; // X
    LocalEdgeVector<2> vel_gauss = evy.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_upper.array() * (vel_gauss.array()).max(0));
    phiup.block<1, 6>(c, 0) -= timemesh.dt / mesh.hy * tmp * BiG23_2;
}

template <int DGdegree>
void boundary_left(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evx,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= timemesh.dt / mesh.hx * std::max(0., -evx(e, 0)) * phi(c, 0);
}
template <>
void boundary_left(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evx,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<1> phi_left
        = LocalEdgeVector<1>(phi(c, 0) - 0.5 * phi(c, 1), phi(c, 2)) * BiGe12;

    // LocalEdgeVector<2> phi_gauss = phi_left * BiGe23; // Y
    LocalEdgeVector<1> vel_gauss = evx.block<1, 2>(e, 0) * BiGe12;
    LocalEdgeVector<1> tmp = (phi_left.array() * (-vel_gauss.array()).max(0));
    phiup.block<1, 3>(c, 0) -= timemesh.dt / mesh.hx * tmp * BiG12_3;
}
template <>
void boundary_left(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evx,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<2> phi_left
        = LocalEdgeVector<2>(phi(c, 0) - 0.5 * phi(c, 1) + 1. / 6. * phi(c, 3),
              phi(c, 2) - 0.5 * phi(c, 5),
              phi(c, 4))
        * BiGe23;
    // LocalEdgeVector<2> phi_gauss = phi_left * BiGe23; // Y
    LocalEdgeVector<2> vel_gauss = evx.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_left.array() * (-vel_gauss.array()).max(0));
    phiup.block<1, 6>(c, 0) -= timemesh.dt / mesh.hx * tmp * BiG23_3;
}

template <int DGdegree>
void boundary_right(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evx,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= timemesh.dt / mesh.hx * std::max(0., evx(e, 0)) * phi(c, 0);
}
template <>
void boundary_right(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evx,
    const size_t c,
    const size_t e)
{
    //    return;
    LocalEdgeVector<1> phi_right
        = LocalEdgeVector<1>(phi(c, 0) + 0.5 * phi(c, 1), phi(c, 2)) * BiGe12;

    // LocalEdgeVector<2> phi_gauss = phi_right * BiGe23; // Y
    LocalEdgeVector<1> vel_gauss = evx.block<1, 2>(e, 0) * BiGe12;
    LocalEdgeVector<1> tmp = (phi_right.array() * (vel_gauss.array().max(0)));

    phiup.block<1, 3>(c, 0) -= timemesh.dt / mesh.hx * tmp * BiG12_1;
}
template <>
void boundary_right(const Mesh& mesh,
    const TimeMesh& timemesh,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evx,
    const size_t c,
    const size_t e)
{
    //    return;
    LocalEdgeVector<2> phi_right // phi in Gauss points on right element
        = LocalEdgeVector<2>(phi(c, 0) + 0.5 * phi(c, 1) + 1. / 6. * phi(c, 3),
              phi(c, 2) + 0.5 * phi(c, 5),
              phi(c, 4))
        * BiGe23;
    // LocalEdgeVector<2> phi_gauss = phi_right * BiGe23; // Y
    LocalEdgeVector<2> vel_gauss = evx.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_right.array() * (vel_gauss.array().max(0)));
    phiup.block<1, 6>(c, 0) -= timemesh.dt / mesh.hx * tmp * BiG23_1;
}

//////////////////////////////////////////////////

//! Computes the "something" like the inverse of the mass matrix. GIVE DETAILS!
template <int DGdegree>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<DGdegree>& inversemasscell);
template <>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<0>& inversemasscell)
{
    inversemasscell(0) = 0; // not used for DG0
}

//  basis functions on (0,h)^2
//
//  1           (0   ,   0)
//  x/h - 1/2   (1/h ,   0)
//  y/h - 1/2   (0   , 1/2)

template <>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<1>& inversemasscell)
{
    // Mass = h^2 diag(1.0, 1.0/12, 1.0/12) / k

    // (Psi[0], nabla Phi) = diag(0, h, h)

    // scaling comes from integral of cell term divided by integral of mass term
    // scaling from derivative is added later in the computation
    inversemasscell(0) = 0.0; // used for cell-term (v Psi, phi)
    inversemasscell(1) = timemesh.dt * 12.0;
    inversemasscell(2) = timemesh.dt * 12.0;
}

template <>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<2>& inversemasscell)
{
    // Mass = h^2 diag(1.0, 1.0/12, 1.0/12) / k

    // (Psi[0], nabla Phi) = diag(0, h, h)

    inversemasscell(0) = 0.0; // used for cell-term (v Psi, phi)
    inversemasscell(1) = timemesh.dt * 12.0;
    inversemasscell(2) = timemesh.dt * 12.0;
    inversemasscell(3) = timemesh.dt * 180.0;
    inversemasscell(4) = timemesh.dt * 180.0;
    inversemasscell(5) = timemesh.dt * 144.0;
}

//! Applies the linear transport operator 'div (v Psi)' in upwind formulation,
//! scaled with inverse mass matrix and with the time step
//! The result is written to phiup (which is set to zero at the beginning)
template <int DGdegree>
void transportoperator(const Mesh& mesh,
    const TimeMesh& timemesh,
    const CellVector<DGdegree>& vx,
    const CellVector<DGdegree>& vy,
    const EdgeVector<DGdegree>& evx,
    const EdgeVector<DGdegree>& evy,
    const CellVector<DGdegree>& phi,
    CellVector<DGdegree>& phiup)
{
    phiup.zero();

    LocalCellVector<DGdegree> inversemasscell;
    inversemassmatrix(mesh, timemesh, inversemasscell);

    GlobalTimer.start("-- -- --> cell term");
    // Cell terms
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx;
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic)
            cell_term<DGdegree>(mesh, inversemasscell, phiup, phi, vx, vy, ic);
    }
    GlobalTimer.stop("-- -- --> cell term");

    GlobalTimer.start("-- -- --> edge terms");
    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row
        size_t ie = iy * (mesh.nx + 1) + 1; // first index of inner velocity in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic, ++ie)
            edge_term_Y<DGdegree>(
                mesh, timemesh, phiup, phi, evx, ic, ic + 1, ie);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        size_t ie = ix + mesh.nx; // first index of inner velocity in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx, ie += mesh.nx)
            edge_term_X<DGdegree>(
                mesh, timemesh, phiup, phi, evy, ic, ic + mesh.nx, ie);
    }
    GlobalTimer.stop("-- -- --> edge terms");

    // boundaries
    GlobalTimer.start("-- -- --> boundaries");
    // lower & upper
    //#pragma omp parallel for
    const size_t eupper0 = mesh.nx * mesh.ny;
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        const size_t clower = ix;
        const size_t elower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;
        const size_t eupper = eupper0 + ix;

        boundary_lower<DGdegree>(
            mesh, timemesh, phiup, phi, evy, clower, elower);
        boundary_upper<DGdegree>(
            mesh, timemesh, phiup, phi, evy, cupper, eupper);
    }
    // left & right
    //#pragma omp parallel for
    size_t eright0 = mesh.nx;

    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t eleft = iy * (mesh.nx + 1);
        const size_t cright = eright0 - 1 + iy * mesh.nx;
        const size_t eright = eright0 + iy * (mesh.nx + 1);

        boundary_left<DGdegree>(
            mesh, timemesh, phiup, phi, evx, cleft, eleft);
        boundary_right<DGdegree>(
            mesh, timemesh, phiup, phi, evx, cright, eright);
    }

    GlobalTimer.stop("-- -- --> boundaries");
}

} // namespace Nextsim

/*----------------------------   dgtimestepping.h ---------------------------*/
/* end of #ifndef __dgtimestepping_H */
#endif
/*----------------------------   dgtimestepping.h ---------------------------*/
