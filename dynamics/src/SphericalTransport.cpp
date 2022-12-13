/*!
 * @file SphericalTransport.cpp
 * @date July 10, 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "SphericalTransport.hpp"
#include "ParametricTools.hpp"
//#include "dgTimeStepping.hpp"
#include "Interpolations.hpp"
#include "stopwatch.hpp"

namespace Nextsim {

extern Timer GlobalTimer;

#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))

//! returns the localization of the cell vector to the edges
template <int DG>
Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)>
leftedgeofcell(const DGVector<DG>& cv, size_t eid);
template <int DG>
Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)>
rightedgeofcell(const DGVector<DG>& cv, size_t eid);
template <int DG>
Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)>
bottomedgeofcell(const DGVector<DG>& cv, size_t eid);
template <int DG>
Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)>
topedgeofcell(const DGVector<DG>& cv, size_t eid);

// dG0 (1 in cell, 1 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
leftedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
rightedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
bottomedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
topedgeofcell(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}

// dG1 (3 in cell, 2 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
leftedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
rightedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
bottomedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
topedgeofcell(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
leftedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5), cv(eid, 4));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
rightedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5), cv(eid, 4));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
bottomedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5), cv(eid, 3));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
topedgeofcell(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5), cv(eid, 3));
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Computes the normal velocity on each edge
  // The velocity is scaled with the length of the edge, therefore
  // the length-element J is already included!
template <int DG>
void SphericalTransport<DG>::reinitnormalvelocity()
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
	  normalvel_Y.row(ey) += 0.5 * (tangent_left(0, 1) * leftedgeofcell<DG>(velx, cy) - tangent_left(0, 0) * leftedgeofcell<DG>(vely, cy));
	  
	  
	  // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
	  const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_right = smesh.edgevector(ey + 1, ey + smesh.nx + 2);
	  normalvel_Y.row(ey + 1) += 0.5 * (tangent_right(0, 1) * rightedgeofcell<DG>(velx, cy) - tangent_right(0, 0) * rightedgeofcell<DG>(vely, cy));
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

            normalvel_X.row(cx) += 0.5 * (-tangent_bottom(0, 1) * bottomedgeofcell<DG>(velx, cx) + tangent_bottom(0, 0) * bottomedgeofcell<DG>(vely, cx));

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_top = smesh.edgevector(nx + smesh.nx + 1, nx + smesh.nx + 2);

            normalvel_X.row(cx + smesh.nx) += 0.5 * (-tangent_top(0, 1) * topedgeofcell<DG>(velx, cx) + tangent_top(0, 0) * topedgeofcell<DG>(vely, cx));
        }

        // scale boundary
        normalvel_X.row(ix) *= 2.0;
        normalvel_X.row(ix + smesh.ny * smesh.nx) *= 2.0;
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
void SphericalTransport<DG>::prepareAdvection(const CGVector<CG>& cg_vx, const CGVector<CG>& cg_vy)
{
    Nextsim::Interpolations::CG2DGSpherical(smesh, GetVx(), cg_vx);
    Nextsim::Interpolations::CG2DGSpherical(smesh, GetVy(), cg_vy);
    reinitnormalvelocity();
}

////////////////////////////////////////////////// CELL TERM

template <int DG>
void cell_term(const ParametricMesh& smesh, double dt,
    DGVector<DG>& phiup, const DGVector<DG>& phi,
    const DGVector<DG>& vx,
    const DGVector<DG>& vy, const size_t ic);
template <>
void cell_term(const ParametricMesh& smesh, double dt,
    DGVector<1>& phiup, const DGVector<1>& phi,
    const DGVector<1>& vx,
    const DGVector<1>& vy, const size_t ic) { }

// higher order terms with gauss quadrature
//
//     // - (A v, \nabla PSI)
//   = - R  ( A  v_lon d_lon psi + v_lat cos(lat) d_lat phi) 
// -  wq *  phi[q] * v[q] * ( Jq * JT^{-T} ) [BIGx, BIGy]

template <int DG>
void cell_term(const ParametricMesh& smesh, double dt,
    DGVector<DG>& phiup,
    const DGVector<DG>& phi,
    const DGVector<DG>& vx,
    const DGVector<DG>& vy, const size_t eid)
{
  if (smesh.landmask[eid]==0)
    return;
  
#define NGP ( (DG <=3 ? 2 : 3 ) )

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
    const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> cos_lat = (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(smesh, eid).row(1).array()*M_PI/180.).cos();


    const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dxT = ParametricTools::dxT<NGP>(smesh, eid);
    const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dyT = ParametricTools::dyT<NGP>(smesh, eid);

    // [J dT^{-T} nabla phi]_1
    phiup.row(eid) += dt * ((PSIx<DG, NGP>.array().rowwise() * dyT.row(1).array() - PSIy<DG, NGP>.array().rowwise() * dxT.row(1).array()).rowwise() * vx_gauss.array() +
			    // [J dT^{-T} nabla phi]_2
			    (PSIy<DG, NGP>.array().rowwise() * dxT.row(0).array() - PSIx<DG, NGP>.array().rowwise() * dyT.row(0).array()).rowwise() * (vy_gauss.array() * cos_lat.array()))
                               .matrix()
        * phi_gauss.transpose();

#undef NGP
}
////////////////////////////////////////////////// BOUNDARY HANDLING


template <int DG>
void boundary_lower(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c, const size_t e)
{
    // GP = DGEDGE
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    // block<1, 2>(e, 0) * PSIe<2,2>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((bottomedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (-vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 0>;
}
template <int DG>
void boundary_upper(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c, const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((topedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 2>;
}
template <int DG>
void boundary_left(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c, const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((leftedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (-vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 3>;
}
template <int DG>
void boundary_right(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c, const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((rightedgeofcell<DG>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (vel_gauss.array().max(0)));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 1>;
}



////////////////////////////////////////////////// EDGE TERMS

inline void edge_term_X(const ParametricMesh& smesh, const double dt, DGVector<1>& phiup, const DGVector<1>& phi, // DG0 (1)
    const EdgeVector<1>& normalvel_X, const size_t c1, const size_t c2, const size_t ie)
{
  if (smesh.landmask[c1]==0) return;
  if (smesh.landmask[c2]==0) return;
  
    double bottom = phi(c1, 0);
    double top = phi(c2, 0);
    double vel = normalvel_X(ie, 0);

    const Eigen::Matrix<Nextsim::FloatType, 2, 2> coords = smesh.coordinatesOfEdgeX(ie);
    const double lat = 0.5*(coords(0,1) + coords(0,1))/180.0*M_PI;
    double Jlon = SQR(coords(0,0)-coords(1,0));
    double Jlat = SQR(coords(0,1)-coords(1,1));
    double JJ   = Jlon + Jlat;
    Jlon /= JJ;
    Jlat /= JJ;
    double ds = sqrt(Jlon * cos(lat)*cos(lat) + Jlat);
    
    phiup(c1, 0) -= dt * ds * (std::max(vel, 0.) * bottom + std::min(vel, 0.) * top);
    phiup(c2, 0) += dt * ds * (std::max(vel, 0.) * bottom + std::min(vel, 0.) * top);
}
inline void edge_term_Y(const ParametricMesh& smesh, const double dt, DGVector<1>& phiup, const DGVector<1>& phi, // DG0 (1)
    const EdgeVector<1>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie)
{
  if (smesh.landmask[c1]==0) return;
  if (smesh.landmask[c2]==0) return;
  
    double left = phi(c1, 0);
    double right = phi(c2, 0);
    double vel = normalvel_Y(ie, 0);

    
    const Eigen::Matrix<Nextsim::FloatType, 2, 2> coords = smesh.coordinatesOfEdgeY(ie);
    const double lat = 0.5*(coords(0,1) + coords(0,1))/180.0*M_PI;
    double Jlon = SQR(coords(0,0)-coords(1,0));
    double Jlat = SQR(coords(0,1)-coords(1,1));
    double JJ   = Jlon + Jlat;
    Jlon /= JJ;
    Jlat /= JJ;
    double ds = sqrt(Jlon * cos(lat)*cos(lat) + Jlat);
 
    phiup(c1, 0) -= dt * ds * (std::max(vel, 0.) * left + std::min(vel, 0.) * right);
    phiup(c2, 0) += dt * ds * (std::max(vel, 0.) * left + std::min(vel, 0.) * right);
}

template <int DG>
inline void edge_term_X(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup, const DGVector<DG>& phi, // DG1 (3)
    const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c1, const size_t c2, const size_t ie)
{
  if (smesh.landmask[c1]==0) return;
  if (smesh.landmask[c2]==0) return;
  
  const LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(ie) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;

  // for integration get the rel. length of the edge in lon/lat directions (squared)
  const Eigen::Matrix<Nextsim::FloatType, 2, 2> coords = smesh.coordinatesOfEdgeX(ie);
  double Jlon = SQR(coords(0,0)-coords(1,0));
  double Jlat = SQR(coords(0,1)-coords(1,1));
  double JJ   = Jlon + Jlat;
  Jlon /= JJ;
  Jlat /= JJ;
 
  const Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)> cos_lat = (ParametricTools::getGaussPointsOnEdgeX<EDGEDOFS(DG)>(smesh,ie).row(1).array()/180.*M_PI).cos();
  const Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)> ds_gauss = (Jlat + Jlon * cos_lat.array().square()).sqrt();

  const LocalEdgeVector<EDGEDOFS(DG)> tmp =//ds_gauss.array() *
    (vel_gauss.array().max(0) * (topedgeofcell<DG>   (phi, c1) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() + 
     vel_gauss.array().min(0) * (bottomedgeofcell<DG>(phi, c2) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() );
    
    phiup.row(c1) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 2>;
    phiup.row(c2) += dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 0>;
}
template <int DG>
inline void edge_term_Y(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup, const DGVector<DG>& phi, // DG1 (3)
    const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie)
{
  if (smesh.landmask[c1]==0) return;
  if (smesh.landmask[c2]==0) return;
  
    const LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(ie) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;

    // for integration get the rel. length of the edge in lon/lat directions (squared)
    const Eigen::Matrix<Nextsim::FloatType, 2, 2> coords = smesh.coordinatesOfEdgeY(ie);
    double Jlon = SQR(coords(0,0)-coords(1,0));
    double Jlat = SQR(coords(0,1)-coords(1,1));
    double JJ   = Jlon + Jlat;
    Jlon /= JJ;
    Jlat /= JJ;

    const Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)> cos_lat = (ParametricTools::getGaussPointsOnEdgeX<EDGEDOFS(DG)>(smesh,ie).row(1).array()/180*M_PI).cos();
    const Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)> ds_gauss = (Jlat + Jlon * cos_lat.array().square()).sqrt();
	
    
    const LocalEdgeVector<EDGEDOFS(DG)> tmp =// ds_gauss.array() *
      (vel_gauss.array().max(0) * (rightedgeofcell<DG>(phi, c1) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() +
       vel_gauss.array().min(0) * (leftedgeofcell <DG>(phi, c2) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array());

    // - [[psi]] sind we're on the left side
    phiup.row(c1) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 1>;
    phiup.row(c2) += dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 3>;
}

template <int DG>
void SphericalTransportOperator(const ParametricMesh& smesh, const double dt,
    const DGVector<DG>& vx,
    const DGVector<DG>& vy,
    const EdgeVector<EDGEDOFS(DG)>& normalvel_X,
    const EdgeVector<EDGEDOFS(DG)>& normalvel_Y,
    const DGVector<DG>& phi, DGVector<DG>& phiup)
{
    phiup.zero();

    GlobalTimer.start("-- -- --> cell term");
    // Cell terms
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        cell_term<DG>(smesh, dt, phiup, phi, vx, vy, eid);
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


    // Periodic
    for (size_t pc = 0 ; pc<smesh.periodic.size(); ++pc)
      {
#pragma omp parallel for
	for (size_t i = 0; i<smesh.periodic[pc].size();++i)
	  {
	    if (smesh.periodic[pc][i][0]==0) // X-edge (bottom / top)
	      edge_term_X(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc][i][1], smesh.periodic[pc][i][2], smesh.periodic[pc][i][3]);
	    else if (smesh.periodic[pc][i][0]==1) // Y-edge (left / right)
	      edge_term_Y(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc][i][1], smesh.periodic[pc][i][2], smesh.periodic[pc][i][3]);
	    else
	      {
		std::cerr << "Wrong periodic boundary information in the mesh. Boundary side " << smesh.periodic[pc][i][0] << " not valid" << std::endl;
		abort();
	      }
	  }
      }
    
    // Dirichlet
    for (size_t seg = 0 ; seg<4; ++seg) // run over the 4 segments (bot, right, top, left)
      {
#pragma omp parallel for
	for (size_t i = 0; i<smesh.dirichlet[seg].size();++i)
	  {
	    const size_t eid = smesh.dirichlet[seg][i];
	    const size_t ix = eid % smesh.nx; // compute 'coordinate' of element
	    const size_t iy = eid / smesh.nx;
    
	    if (seg==0) // bottom
	      boundary_lower(smesh, dt, phiup, phi, normalvel_X, eid, smesh.nx*iy+ix);
	    else if (seg==1) // right
	      boundary_right(smesh, dt, phiup, phi, normalvel_Y, eid, (smesh.nx+1)*iy+ix+1);
	    else if (seg==2) // top
	      boundary_upper(smesh, dt, phiup, phi, normalvel_X, eid, smesh.nx*(iy+1)+ix);
	    else if (seg==3) // left
	      boundary_left (smesh, dt, phiup, phi, normalvel_Y, eid, (smesh.nx+1)*iy+ix);
	    else
	      {
		std::cerr << "Wrong Dirichlet boundary information in the mesh. Boundary side " << smesh.dirichlet[seg][i] << " not valid" << std::endl;
		abort();
	      }
	  }
      }
    
    GlobalTimer.stop("-- -- --> boundaries");

    ////// APPLY INVERSE MASS MATRIX !!! Switch to precomputed!
    
    GlobalTimer.start("-- -- --> inverse mass");
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        phiup.row(eid) =  (SphericalTools::massMatrix<DG>(smesh, eid).inverse() * phiup.row(eid).transpose());
    }

    GlobalTimer.stop("-- -- --> inverse mass");
}

template <int DG>
void SphericalTransport<DG>::step_rk1(const double dt, DGVector<DG>& phi)
{
    SphericalTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1);

    phi += tmp1;
}

template <int DG>
void SphericalTransport<DG>::step_rk2(const double dt, DGVector<DG>& phi)
{
    SphericalTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1); // tmp1 = k * F(u)

    phi += tmp1; // phi = phi + k * F(u)     (i.e.: implicit Euler)

    SphericalTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp2); // tmp1 = k * F( u + k * F(u) )

    phi += 0.5 * (tmp2 - tmp1);
}

template <int DG>
void SphericalTransport<DG>::step_rk3(const double dt, DGVector<DG>& phi)
{
    SphericalTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
        tmp1); // tmp1 = k * F(u)  // K1 in Heun(3)
    tmp1 += phi; // phi + h f(phi)

    SphericalTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, tmp1,
        tmp2); // tmp2 = f( u + h f(u) )
    tmp2 += tmp1;
    tmp2 *= 0.25;
    tmp2 += 0.75 * phi;

    SphericalTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, tmp2,
        tmp3); // k * F(k1) // K2 in Heun(3)
    tmp3 += tmp2;

    phi *= 1.0 / 3.0;
    phi += 2.0 / 3.0 * tmp3;

    // parametricTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
    //     tmp1); // tmp1 = k * F(u)  // K1 in Heun(3)

    // phi += 1. / 3. * tmp1; // phi = phi + k/3 * F(u)   (i.e.: implicit Euler)
    // parametricTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
    //     tmp2); // k * F(k1) // K2 in Heun(3)
    // phi -= 1. / 3. * tmp1; // phi = phi + k/3 * F(u)   (i.e.: implicit Euler)

    // phi += 2. / 3. * tmp2;
    // parametricTransportOperator<DG>(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi,
    //     tmp3); // k * F(k2) // K3 in Heun(3)
    // phi -= 2. / 3. * tmp2;

    // phi += 0.25 * tmp1 + 0.75 * tmp3;
}

template <int DG>
void SphericalTransport<DG>::step(const double dt, DGVector<DG>& phi)
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

#undef EDGEDOFS

template class SphericalTransport<1>;
template class SphericalTransport<3>;
template class SphericalTransport<6>;

template void SphericalTransport<1>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void SphericalTransport<1>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void SphericalTransport<3>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void SphericalTransport<3>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void SphericalTransport<6>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void SphericalTransport<6>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);

} /* namespace Nextsim */
