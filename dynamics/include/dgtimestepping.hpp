/*----------------------------   dgtimestepping.h ---------------------------*/
#ifndef __dgtimestepping_H
#define __dgtimestepping_H
/*----------------------------   dgtimestepping.h ---------------------------*/

#include "dgvector_manipulations.hpp"
#include "timemesh.hpp"

namespace Nextsim {

// Zero's for proper integration
//
// dg0 - none
//
// dg1 - a0 + a1*(x-1/2)
// a1 != 0  ==>  z = 1/2 - a0/a1
//
// dg2 - a0 + a1*(x-1/2) + a2 * ( (x-1/2)^2 - 1/12 )
//
// disc:
// a2 !=0 sonst dg1
// D = -36 * a0 * a2 + 9 * a1 * a1 + 3 * a2 * a2
// D>0 ==> z0/1 = ( a2 - a1 +/- sqrt(D) ) / (2 * a2)
//

// pre-computed matrices for assembling dG-transport
// the Gauss rule for dG(n) is set to n+1 points
// this might not be enough?
#include "dgbasisfunctions_gausspoints.hpp"

// compute the cell terms starting dG(1)
template <int DGdegree>
void cell_term(const LocalCellVector<DGdegree> inversemasscell,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const CellVector<DGdegree>& vx,
    const CellVector<DGdegree>& vy,
    const size_t ic);
template <>
void cell_term(const LocalCellVector<0> inversemasscell,
    CellVector<0>& phiup,
    const CellVector<0>& phi,
    const CellVector<0>& vx,
    const CellVector<0>& vy,
    const size_t ic)
{
}
template <>
void cell_term(const LocalCellVector<1> inversemasscell,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const CellVector<1>& vx,
    const CellVector<1>& vy,
    const size_t ic)
{
    // - (Psi v, \nabla phi)

    // vx, vy dG(0) !!!
    // nable phi = (0, 1, 1)

    // AV[0] * dx Phi(1) * VX[0]
    phiup(ic, 1) += inversemasscell(1) * (phi(ic, 0) * vx(ic, 0) + 1. / 12. * phi(ic, 1) * vx(ic, 1) + 1. / 12. * phi(ic, 2) * vx(ic, 2));
    // AV[0] * dy Phi(2) * VX[0]
    phiup(ic, 2) += inversemasscell(2) * (phi(ic, 0) * vy(ic, 0) + 1. / 12. * phi(ic, 1) * vy(ic, 1) + 1. / 12. * phi(ic, 2) * vy(ic, 2));

    // alle anderen Beitraege 0
}
template <>
void cell_term(const LocalCellVector<2> inversemasscell,
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

    // AV[0] * dx Phi(1) * VX[0]

    phiup(ic, 1) += inversemasscell(1) * (phi(ic, 0) * vx(ic, 0) + 1. / 12. * phi(ic, 1) * vx(ic, 1) + 1. / 12. * phi(ic, 2) * vx(ic, 2) + 1. / 180. * phi(ic, 3) * vx(ic, 3) + 1. / 180. * phi(ic, 4) * vx(ic, 4) + 1. / 144. * phi(ic, 5) * vx(ic, 5));

    phiup(ic, 2) += inversemasscell(2) * (phi(ic, 0) * vy(ic, 0) + 1. / 12. * phi(ic, 1) * vy(ic, 1) + 1. / 12. * phi(ic, 2) * vy(ic, 2) + 1. / 180. * phi(ic, 3) * vy(ic, 3) + 1. / 180. * phi(ic, 4) * vy(ic, 4) + 1. / 144. * phi(ic, 5) * vy(ic, 5));

    phiup(ic, 3) += inversemasscell(3) * (1. / 6. * (phi(ic, 1) * vx(ic, 0) + phi(ic, 0) * vx(ic, 1)) + 1. / 90. * (phi(ic, 3) * vx(ic, 1) + phi(ic, 1) * vx(ic, 3)) + 1. / 72. * (phi(ic, 2) * vx(ic, 5) + phi(ic, 5) * vx(ic, 2)));

    phiup(ic, 4) += inversemasscell(4) * (1. / 6. * (phi(ic, 2) * vy(ic, 0) + phi(ic, 0) * vy(ic, 2)) + 1. / 90. * (phi(ic, 2) * vy(ic, 4) + phi(ic, 4) * vy(ic, 2)) + 1. / 72. * (phi(ic, 1) * vy(ic, 5) + phi(ic, 5) * vy(ic, 1)));

    phiup(ic, 5) += inversemasscell(5) * (1. / 12. * (phi(ic, 0) * vx(ic, 2) + phi(ic, 2) * vx(ic, 0)) + 1. / 144. * (phi(ic, 1) * vx(ic, 5) + phi(ic, 5) * vx(ic, 1)) + 1. / 180. * (phi(ic, 2) * vx(ic, 4) + phi(ic, 4) * vx(ic, 2)));

    phiup(ic, 5) += inversemasscell(5) * (1. / 12. * (phi(ic, 0) * vy(ic, 1) + phi(ic, 1) * vy(ic, 0)) + 1. / 144. * (phi(ic, 2) * vy(ic, 5) + phi(ic, 5) * vy(ic, 2)) + 1. / 180. * (phi(ic, 1) * vy(ic, 3) + phi(ic, 3) * vy(ic, 1)));
    // alle anderen Beitraege 0
}

// compute the edge terms for the horizontal edges:  n = (0, +/- 1)
template <int DGdegree>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<DGdegree> scale,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evy,
    const size_t c1,
    const size_t c2,
    const size_t e);

template <>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<0> scale,
    CellVector<0>& phiup,
    const CellVector<0>& phi,
    const EdgeVector<0>& evy,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    const double avgphi = 0.5 * (phi(c1, 0) + phi(c2, 0));
    const double jmpphi = phi(c1, 0) - phi(c2, 0);

    // + < {{Psi}, v* [[phi]] > = h * {{Psi}} * (vright - vleft)
    const double pa = scale(0, 0) * avgphi * evy(ie, 0);
    phiup(c1, 0) -= pa;
    phiup(c2, 0) += pa;

    // + 1/2 < |v n|  [[Psi]] [[phi]] = h/2 * |v n| [[Psi]]  (right - left)
    const double pj = 0.5 * scale(0, 0) * jmpphi * fabs(evy(ie, 0));
    phiup(c1, 0) -= pj;
    phiup(c2, 0) += pj;
}

template <>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<1> scale,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evy,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    // average. cell: (1, x/h-1/2, y/h-1/2)  edge (1, x/h-1/2)
    LocalEdgeVector<1> avgphi(0.5 * (phi(c1, 0) + phi(c2, 0)) + 0.25 * (phi(c1, 2) - phi(c2, 2)),
        0.5 * (phi(c1, 1) + phi(c2, 1)));
    LocalEdgeVector<1> jmpphi(phi(c1, 0) + 0.5 * phi(c1, 2) - (phi(c2, 0) - 0.5 * phi(c2, 2)),
        phi(c1, 1) - phi(c2, 1));

    LocalEdgeVector<1> avg_gauss = avgphi * BiGe12;

    LocalEdgeVector<1> jmp_gauss = jmpphi * BiGe12;
    LocalEdgeVector<1> vel_gauss = evy.block<1, 2>(ie, 0) * BiGe12;
    LocalEdgeVector<1> tmp = (avg_gauss.array() * vel_gauss.array() + 0.5 * jmp_gauss.array() * vel_gauss.array().abs());

    phiup.block<1, 3>(c1, 0) -= timemesh.k / mesh.h * tmp * BiG12_2;
    phiup.block<1, 3>(c2, 0) += timemesh.k / mesh.h * tmp * BiG12_0;
}

template <>
void edge_term_X(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<2> scale,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evy,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    //             2
    // 1  Y  2     X
    //             1
    //
    //
    // X1
    //                     1     x-1/2     (x-1/2)^2-1/12
    //
    // 1                   1     0
    // x-1/2               0     1
    // y-1/2               1/2   0
    // (x-1/2)^2-1/12      0     0         1
    // (y-1/2)^2-1/12      1/6   0
    // (x-1/2)(y-1/2)      0     1/2

    // X2
    //                     1     x-1/2     (x-1/2)^2-1/12
    //
    // 1                   1     0
    // x-1/2               0     1
    // y-1/2              -1/2   0
    // (x-1/2)^2-1/12      0     0         1
    // (y-1/2)^2-1/12      1/6   0
    // (x-1/2)(x-1/12)     0    -1/2

    // average. cell: (1, x/h-1/2, y/h-1/2)  edge (1, x/h-1/2)
    LocalEdgeVector<2> avgphi(
        0.5 * ((phi(c1, 0) + 0.5 * phi(c1, 2) + 1. / 6. * phi(c1, 4)) + (phi(c2, 0) - 0.5 * phi(c2, 2) + 1. / 6. * phi(c2, 4))),
        0.5 * ((phi(c1, 1) + 0.5 * phi(c1, 5)) + (phi(c2, 1) - 0.5 * phi(c2, 5))),
        0.5 * ((phi(c1, 3)) + (phi(c2, 3))));
    LocalEdgeVector<2> jmpphi(
        ((phi(c1, 0) + 0.5 * phi(c1, 2) + 1. / 6. * phi(c1, 4)) - (phi(c2, 0) - 0.5 * phi(c2, 2) + 1. / 6. * phi(c2, 4))),
        ((phi(c1, 1) + 0.5 * phi(c1, 5)) - (phi(c2, 1) - 0.5 * phi(c2, 5))),
        ((phi(c1, 3)) - (phi(c2, 3))));

    LocalEdgeVector<2> avg_gauss = avgphi * BiGe23;
    LocalEdgeVector<2> jmp_gauss = jmpphi * BiGe23;
    LocalEdgeVector<2> vel_gauss = evy.block<1, 3>(ie, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (avg_gauss.array() * vel_gauss.array() + 0.5 * jmp_gauss.array() * vel_gauss.array().abs());

    phiup.block<1, 6>(c1, 0) -= timemesh.k / mesh.h * tmp * BiG23_2;
    phiup.block<1, 6>(c2, 0) += timemesh.k / mesh.h * tmp * BiG23_0;
}

// compute the edge terms for the vertical edges:  n = (+/- 1, 0)
template <int DGdegree>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<DGdegree> scale,
    CellVector<DGdegree>& phiup,

    const CellVector<DGdegree>& phi,

    const EdgeVector<DGdegree>& evx,
    const size_t c1,
    const size_t c2,
    const size_t e);

template <>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<0> scale,
    CellVector<0>& phiup,

    const CellVector<0>& phi,

    const EdgeVector<0>& evx,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    const double avgphi = 0.5 * (phi(c1, 0) + phi(c2, 0));
    const double jmpphi = phi(c1, 0) - phi(c2, 0);

    // + < {{Psi}, v* [[phi]] > = h * {{Psi}} * (vright - vleft)
    const double pa = scale(0) * avgphi * evx(ie, 0);
    phiup(c1, 0) -= pa;
    phiup(c2, 0) += pa;

    // + 1/2 < |v n|  [[Psi]] [[phi]] = h/2 * |v n| [[Psi]]  (right - left)
    const double pj = 0.5 * scale(0) * jmpphi * fabs(evx(ie, 0));
    phiup(c1, 0) -= pj;
    phiup(c2, 0) += pj;
}

template <>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<1> scale,
    CellVector<1>& phiup,
    const CellVector<1>& phi,
    const EdgeVector<1>& evx,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{
    // average. cell: (1, x-1/2, y-1/2)  edge (1, y-1/2)
    LocalEdgeVector<1> avgphi(0.5 * (phi(c1, 0) + phi(c2, 0)) + 0.25 * (phi(c1, 1) - phi(c2, 1)),
        0.5 * (phi(c1, 2) + phi(c2, 2)));
    LocalEdgeVector<1> jmpphi(phi(c1, 0) + 0.5 * phi(c1, 1) - (phi(c2, 0) - 0.5 * phi(c2, 1)),
        phi(c1, 2) - phi(c2, 2));

    LocalEdgeVector<1> avg_gauss = avgphi * BiGe12;
    LocalEdgeVector<1> jmp_gauss = jmpphi * BiGe12;
    LocalEdgeVector<1> vel_gauss = evx.block<1, 2>(ie, 0) * BiGe12;
    LocalEdgeVector<1> tmp = (avg_gauss.array() * vel_gauss.array() + 0.5 * jmp_gauss.array() * vel_gauss.array().abs());

    // - [[psi]] sind we're on the left side
    phiup.block<1, 3>(c1, 0) -= timemesh.k / mesh.h * tmp * BiG12_1;
    phiup.block<1, 3>(c2, 0) += timemesh.k / mesh.h * tmp * BiG12_3;
}

template <>
void edge_term_Y(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<2> scale,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evx,
    const size_t c1,
    const size_t c2,
    const size_t ie)
{

    //             2
    // 1  Y  2     X
    //             1
    //
    //
    // Y1
    //                     1     y-1/2     (y-1/2)^2-1/12
    //
    // 1                   1     0
    // x-1/2               1/2   0
    // y-1/2               0     1
    // (x-1/2)^2-1/12      1/6   0
    // (y-1/2)^2-1/12      0     0         1
    // (x-1/2)(y-1/2)      0     1/2

    // Y2
    //                     1     y-1/2     (y-1/2)^2-1/12
    //
    // 1                   1     0
    // x-1/2              -1/2   0
    // y-1/2               0     1
    // (x-1/2)^2-1/12      1/6   0
    // (y-1/2)^2-1/12      0     0         1
    // (x-1/2)(y-1/12)     0    -1/2

    // average. cell: (1, x/h-1/2, y/h-1/2)  edge (1, x/h-1/2)
    LocalEdgeVector<2> avgphi(
        0.5 * ((phi(c1, 0) + 0.5 * phi(c1, 1) + 1. / 6. * phi(c1, 3)) + (phi(c2, 0) - 0.5 * phi(c2, 1) + 1. / 6. * phi(c2, 3))),
        0.5 * ((phi(c1, 2) + 0.5 * phi(c1, 5)) + (phi(c2, 2) - 0.5 * phi(c2, 5))),
        0.5 * ((phi(c1, 4)) + (phi(c2, 4))));
    LocalEdgeVector<2> jmpphi(
        ((phi(c1, 0) + 0.5 * phi(c1, 1) + 1. / 6. * phi(c1, 3)) - (phi(c2, 0) - 0.5 * phi(c2, 1) + 1. / 6. * phi(c2, 3))),
        ((phi(c1, 2) + 0.5 * phi(c1, 5)) - (phi(c2, 2) - 0.5 * phi(c2, 5))),
        ((phi(c1, 4)) - (phi(c2, 4))));

    LocalEdgeVector<2> avg_gauss = avgphi * BiGe23;
    LocalEdgeVector<2> jmp_gauss = jmpphi * BiGe23;
    LocalEdgeVector<2> vel_gauss = evx.block<1, 3>(ie, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (avg_gauss.array() * vel_gauss.array() + 0.5 * jmp_gauss.array() * vel_gauss.array().abs());

    // - [[psi]] sind we're on the left side
    phiup.block<1, 6>(c1, 0) -= timemesh.k / mesh.h * tmp * BiG23_1;
    phiup.block<1, 6>(c2, 0) += timemesh.k / mesh.h * tmp * BiG23_3;
}

//////////////////////////////////////////////////

//  <<  <v*n>^+ Psi Phi >>

template <int DGdegree>
void boundary_lower(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<DGdegree> scale,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evy,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= std::max(0., -evy(e, 0)) * phi(c, 0);
}
template <>
void boundary_lower(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<2> scale,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evy,
    const size_t c,
    const size_t e)
{
    // 1, x-1/2, y-1/2, (x-1/2)^2-1/12, (y-1/2)^2-1/12, (x-1/2)(y-1/2)
    LocalEdgeVector<2> phi_lower(phi(c, 0) - 0.5 * phi(c, 2) + 1. / 6. * phi(c, 4),
        phi(c, 1) - 0.5 * phi(c, 5),
        phi(c, 3));
    // LocalEdgeVector<2> phi_gauss = phi_lower * BiGe23; // X
    LocalEdgeVector<2> vel_gauss = evy.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_lower.array() * (-vel_gauss.array()).max(0));
    phiup.block<1, 6>(c, 0) -= timemesh.k / mesh.h * tmp * BiG23_0;
}
template <int DGdegree>
void boundary_upper(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<DGdegree> scale,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evy,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= std::max(0., evy(e, 0)) * phi(c, 0);
}
template <>
void boundary_upper(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<2> scale,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evy,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<2> phi_upper(phi(c, 0) + 0.5 * phi(c, 2) + 1. / 6. * phi(c, 4),
        phi(c, 1) + 0.5 * phi(c, 5),
        phi(c, 3));
    // LocalEdgeVector<2> phi_gauss = phi_upper * BiGe23; // X
    LocalEdgeVector<2> vel_gauss = evy.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_upper.array() * (vel_gauss.array()).max(0));
    phiup.block<1, 6>(c, 0) -= timemesh.k / mesh.h * tmp * BiG23_2;
}

template <int DGdegree>
void boundary_left(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<DGdegree> scale,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evx,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= std::max(0., -evx(e, 0)) * phi(c, 0);
}
template <>
void boundary_left(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<2> scale,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evx,
    const size_t c,
    const size_t e)
{
    LocalEdgeVector<2> phi_left(phi(c, 0) - 0.5 * phi(c, 1) + 1. / 6. * phi(c, 3),
        phi(c, 2) - 0.5 * phi(c, 5),
        phi(c, 4));
    // LocalEdgeVector<2> phi_gauss = phi_left * BiGe23; // Y
    LocalEdgeVector<2> vel_gauss = evx.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_left.array() * (-vel_gauss.array()).max(0));
    phiup.block<1, 6>(c, 0) -= timemesh.k / mesh.h * tmp * BiG23_3;
}

template <int DGdegree>
void boundary_right(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<DGdegree> scale,
    CellVector<DGdegree>& phiup,
    const CellVector<DGdegree>& phi,
    const EdgeVector<DGdegree>& evx,
    const size_t c,
    const size_t e)
{
    phiup(c, 0) -= std::max(0., evx(e, 0)) * phi(c, 0);
}
template <>
void boundary_right(const Mesh& mesh,
    const TimeMesh& timemesh,
    const LocalCellVector<2> scale,
    CellVector<2>& phiup,
    const CellVector<2>& phi,
    const EdgeVector<2>& evx,
    const size_t c,
    const size_t e)
{
    //    return;
    LocalEdgeVector<2> phi_right(phi(c, 0) + 0.5 * phi(c, 1) + 1. / 6. * phi(c, 3),
        phi(c, 2) + 0.5 * phi(c, 5),
        phi(c, 4));
    // LocalEdgeVector<2> phi_gauss = phi_right * BiGe23; // Y
    LocalEdgeVector<2> vel_gauss = evx.block<1, 3>(e, 0) * BiGe23;
    LocalEdgeVector<2> tmp = (phi_right.array() * (vel_gauss.array().max(0)));
    phiup.block<1, 6>(c, 0) -= -timemesh.k / mesh.h * tmp * BiG23_1;
}

//////////////////////////////////////////////////

template <int DGdegree>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<DGdegree>& inversemasscell,
    LocalCellVector<DGdegree>& inversemassedge);
template <>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<0>& inversemasscell,
    LocalCellVector<0>& inversemassedge)
{
    inversemasscell(0) = timemesh.k;

    inversemassedge(0) = timemesh.k / mesh.h;
}

//  basis functions on (0,h)^2
//
//  1           (0   ,   0)
//  x/h - 1/2   (1/h ,   0)
//  y/h - 1/2   (0   , 1/2)

template <>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<1>& inversemasscell,
    LocalCellVector<1>& inversemassedge)
{
    // Mass = h^2 diag(1.0, 1.0/12, 1.0/12) / k

    // (Psi[0], nabla Phi) = diag(0, h, h)

    inversemasscell(0) = 0.0; // used for cell-term (v Psi, phi)
    inversemasscell(1) = timemesh.k / mesh.h * 12.0;
    inversemasscell(2) = timemesh.k / mesh.h * 12.0;

    // MassEdge = diag(h, h/12)

    // Mass^-1 * MassEdge                           // used in edge terms like
    // <v {Psi}, [phi]>
    inversemassedge(0) = timemesh.k / mesh.h; // Psi[0] phi[0]    +/- 0.5 * h for phi[1/2]
    inversemassedge(1) = timemesh.k / mesh.h * 12.0;
    inversemassedge(2) = timemesh.k / mesh.h * 12.0;
}

template <>
void inversemassmatrix(const Mesh& mesh,
    const TimeMesh& timemesh,
    LocalCellVector<2>& inversemasscell,
    LocalCellVector<2>& inversemassedge)
{
    // Mass = h^2 diag(1.0, 1.0/12, 1.0/12) / k

    // (Psi[0], nabla Phi) = diag(0, h, h)

    inversemasscell(0) = 0.0; // used for cell-term (v Psi, phi)
    inversemasscell(1) = timemesh.k / mesh.h * 12.0;
    inversemasscell(2) = timemesh.k / mesh.h * 12.0;
    inversemasscell(3) = timemesh.k / mesh.h * 180.0;
    inversemasscell(4) = timemesh.k / mesh.h * 180.0;
    inversemasscell(5) = timemesh.k / mesh.h * 144.0;

    // MassEdge = diag(h, h/12)

    // Mass^-1 * MassEdge                           // used in edge terms like
    // <v {Psi}, [phi]>
    inversemassedge(0) = timemesh.k / mesh.h; // Psi[0] phi[0]    +/- 0.5 * h for phi[1/2]
    inversemassedge(1) = timemesh.k / mesh.h * 12.0;
    inversemassedge(2) = timemesh.k / mesh.h * 12.0;
    inversemassedge(3) = timemesh.k / mesh.h * 180.0;
    inversemassedge(4) = timemesh.k / mesh.h * 180.0;
    inversemassedge(5) = timemesh.k / mesh.h * 144.0;
}

template <int DGdegree>
void fwdeuler_step(const Mesh& mesh,
    const TimeMesh& timemesh,
    const CellVector<DGdegree>& vx,
    const CellVector<DGdegree>& vy,
    const EdgeVector<DGdegree>& evx,
    const EdgeVector<DGdegree>& evy,
    CellVector<DGdegree>& phi,
    CellVector<DGdegree>& phiup)
{
    phiup.zero();

    LocalCellVector<DGdegree> inversemasscell, inversemassedge;
    inversemassmatrix(mesh, timemesh, inversemasscell, inversemassedge);

    // Cell terms
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx;
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic)
            cell_term<DGdegree>(inversemasscell, phiup, phi, vx, vy, ic);
    }

    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row
        size_t ie = iy * (mesh.nx + 1) + 1; // first index of inner velocity in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic, ++ie)
            edge_term_Y<DGdegree>(
                mesh, timemesh, inversemassedge, phiup, phi, evx, ic, ic + 1, ie);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        size_t ie = ix + mesh.nx; // first index of inner velocity in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx, ie += mesh.nx)
            edge_term_X<DGdegree>(
                mesh, timemesh, inversemassedge, phiup, phi, evy, ic, ic + mesh.nx, ie);
    }

    // boundaries
    // lower & upper
    //#pragma omp parallel for
    const size_t eupper0 = mesh.nx * mesh.ny;
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        const size_t clower = ix;
        const size_t elower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;
        const size_t eupper = eupper0 + ix;

        boundary_lower<DGdegree>(
            mesh, timemesh, inversemassedge, phiup, phi, evy, clower, elower);
        boundary_upper<DGdegree>(
            mesh, timemesh, inversemassedge, phiup, phi, evy, cupper, eupper);
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
            mesh, timemesh, inversemassedge, phiup, phi, evx, cleft, eleft);
        boundary_right<DGdegree>(
            mesh, timemesh, inversemassedge, phiup, phi, evx, cright, eright);
    }
}

} // namespace Nextsim

/*----------------------------   dgtimestepping.h ---------------------------*/
/* end of #ifndef __dgtimestepping_H */
#endif
/*----------------------------   dgtimestepping.h ---------------------------*/
