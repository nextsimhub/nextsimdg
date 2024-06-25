/*!
 * @file DGTransport.cpp
 * @date July 10, 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "DGTransport.hpp"
#include "Interpolations.hpp"
#include "codeGenerationDGinGauss.hpp"

namespace Nextsim {
  
#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))




  /*! 
   * returns the localization of the cell vector to the edges. 
   * EDGE=0 bottom, EDGE=1 right, EDGE=2 top, EDGE=3 left
   */
template <int DG, int EDGE>
Eigen::Matrix<Nextsim::FloatType, 1, EDGEDOFS(DG)>
valueonedge(const DGVector<DG>& cv, size_t eid);

// dG0 (1 in cell, 1 on edge)
template <> // left
Eigen::Matrix<Nextsim::FloatType, 1, 1>
valueonedge<1,3>(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
valueonedge<1,1>(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
valueonedge<1,0>(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 1>
valueonedge<1,2>(const DGVector<1>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 1>(cv(eid, 0));
}

// dG1 (3 in cell, 2 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
valueonedge<3,3>(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 1), cv(eid, 2));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
valueonedge<3,1>(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 1), cv(eid, 2));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
valueonedge<3,0>(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) - 0.5 * cv(eid, 2), cv(eid, 1));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 2>
valueonedge<3,2>(const DGVector<3>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 2>(cv(eid, 0) + 0.5 * cv(eid, 2), cv(eid, 1));
}

// dG2 (6 in cell, 3 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<6,3>(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) - 0.5 * cv(eid, 5), cv(eid, 4));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<6,1>(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
        cv(eid, 2) + 0.5 * cv(eid, 5), cv(eid, 4));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<6,0>(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) - 0.5 * cv(eid, 5), cv(eid, 3));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<6,2>(const DGVector<6>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
        cv(eid, 1) + 0.5 * cv(eid, 5), cv(eid, 3));
}
// dG2+ (8 in cell, 3 on edge)
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<8,3>(const DGVector<8>& cv, size_t eid)
{
    return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
						   cv(eid, 2) - 0.5 * cv(eid, 5) + 1./6. * cv(eid,6),
						   cv(eid, 4) - 0.5 * cv(eid,7));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<8,1>(const DGVector<8>& cv, size_t eid)
{
  return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 1) + 1. / 6. * cv(eid, 3),
						 cv(eid, 2) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 6),
						 cv(eid, 4) + 0.5 * cv(eid, 7));
}
template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<8,0>(const DGVector<8>& cv, size_t eid)
{
  return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) - 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
						 cv(eid, 1) - 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7), 
						 cv(eid, 3) - 0.5 * cv(eid, 6));
}
  template <>
Eigen::Matrix<Nextsim::FloatType, 1, 3>
valueonedge<8,2>(const DGVector<8>& cv, size_t eid) 
{
  return Eigen::Matrix<Nextsim::FloatType, 1, 3>(cv(eid, 0) + 0.5 * cv(eid, 2) + 1. / 6. * cv(eid, 4),
						 cv(eid, 1) + 0.5 * cv(eid, 5) + 1. / 6. * cv(eid, 7),
						 cv(eid, 3) + 0.5 * cv(eid, 6));
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Computes the normal velocity on each edge
  // The velocity is scaled with the length of the edge, therefore
  // the length-element J is already included!
  //
  // The way to compute the normal vector does not differ between
  // the Cartesian and the Sperical version. This is, since the normal
  // is not the actual normal but a scaled version of it. The weight
  // that is missing is corrected in the function edge_term
  // Here, we omit the scaling with
  // cos(theta) and simple use ds = |J| as length scale
  // (instead of ds = (d_theta^2 + cos(theta)^2 d_phi^2)^1/2
  //
template <int DG>
void DGTransport<DG>::reinitnormalvelocity()
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
	  if (smesh.landmask[cy]==0) // skip all land elements
	    continue;
	  
	  // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
	  const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_left = smesh.edgevector(ey, ey + smesh.nx + 1);
	  normalvel_Y.row(ey) += 0.5 * (tangent_left(0, 1) * valueonedge<DG,3>(velx, cy) - tangent_left(0, 0) * valueonedge<DG,3>(vely, cy));
	  
	  // un-normed tangent vector of left edge (pointing up). normal is (y,-x)
	  const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_right = smesh.edgevector(ey + 1, ey + smesh.nx + 2);
	  normalvel_Y.row(ey + 1) += 0.5 * (tangent_right(0, 1) * valueonedge<DG,1>(velx, cy) - tangent_right(0, 0) * valueonedge<DG,1>(vely, cy));
        }

        // scale normal at outer boundaries. Should this also be done at land-masks? Or not relevant sionce always Dirichlet?
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
	  if (smesh.landmask[cx]==0) // skip land elements
	    continue;
	  
            // un-normed tangent vector of bottom edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_bottom = smesh.edgevector(nx, nx + 1);
            normalvel_X.row(cx) += 0.5 * (-tangent_bottom(0, 1) * valueonedge<DG,0>(velx, cx) + tangent_bottom(0, 0) * valueonedge<DG,0>(vely, cx));

            // un-normed tangent vector of top edge (pointing right). normal is (-y,x)
            const Eigen::Matrix<Nextsim::FloatType, 1, 2> tangent_top = smesh.edgevector(nx + smesh.nx + 1, nx + smesh.nx + 2);
            normalvel_X.row(cx + smesh.nx) += 0.5 * (-tangent_top(0, 1) * valueonedge<DG,2>(velx, cx) + tangent_top(0, 0) * valueonedge<DG,2>(vely, cx));
        }

	
	// scale boundary
	normalvel_X.row(ix) *= 2.0;
	normalvel_X.row(ix + smesh.ny * smesh.nx) *= 2.0;
	
    }


    // Not making so much sense... On dirichlet boundaries the velocity is zero.
    
    
//     for (size_t seg = 0 ; seg<4; ++seg) // run over the 4 segments (bot, right, top, left)
//       {
// #pragma omp parallel for
// 	for (size_t i = 0; i<smesh.dirichlet[seg].size();++i)
// 	  {
// 	    const size_t eid = smesh.dirichlet[seg][i];  //! The i of the boundary element 
// 	    const size_t ix = eid % smesh.nx;            //! x & y indices of the element
//            const size_t iy = eid / smesh.nx;
	   
//            if (seg==0) // bottom
//              normalvel_X.row(smesh.nx*iy + ix) *= 2.0;
//            else if (seg==1) // right
//              normalvel_Y.row((smesh.nx+1)*iy + ix+1) *= 2.0;
//            else if (seg==2) // top
//              normalvel_X.row(smesh.nx*(iy+1) + ix) *= 2.0;
//            else if (seg==3) // left
//              normalvel_Y.row((smesh.nx+1)*iy + ix) *= 2.0;
//          }
//       }
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



////////////////////////////////////////////////// CELL TERM

template <>
void DGTransport<1>::cell_term(const ParametricMesh& smesh, double dt,
    DGVector<1>& phiup, const DGVector<1>& phi,
    const DGVector<1>& vx,
    const DGVector<1>& vy, const size_t ic) { }

template <int DG>
void DGTransport<DG>::cell_term(const ParametricMesh& smesh, double dt,
    DGVector<DG>& phiup,
    const DGVector<DG>& phi,
    const DGVector<DG>& vx,
    const DGVector<DG>& vy, const size_t eid)
{
  if (smesh.landmask[eid]==0){
    // TODO remove comments that are here to validate correctness of landmask
    //std::cout << "element id " << eid  << " lm " << smesh.landmask[eid] << " " << phi.row(eid) << std::endl;
    return;
  } else {
    //std::cout << "element id " << eid  << " lm " << smesh.landmask[eid] << " " << phi.row(eid) << std::endl;
  }




  // here, veloicity and phi are evaluated in Gauss points. 
  // vx.row(eid) is a 1 x DG matrix (row-vector), PSI is a matrix of size DG x NQ (no. quad points)
  // * is the Matrix-Matrix Multiplication, result is a 1 x NQ matrix (row vector)
  const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(DG)> vx_gauss = vx.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>; //!< velocity in GP
  const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(DG)> vy_gauss = vy.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>;
  const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(DG)> phi_gauss = phi.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>;

  
  // Das parammap.AdvectionCellTermX/Y[eid] ist die Matrix, die ich Dir aufgeschrieben habe
  // Groesse der Matrix ist jeweils  DG x NQ
  // Das .array() sorgt dafuer, dass die Matrizen nicht mehr als Matrizen aufgefasst werden sondern einfach als Felder. Das Produkt '*' hat
  // dann eine andere Bedeutung
  // Wenn X ein DG x NQ array  ist und Y ein 1 x NQ array, dann bedeutet:
  //    X.rowwise() * Y
  // dass eine Matrix der Groesse DG x NG rauskommt mit (X.rowwise() * Y)_ij = X_ij * Y_j
  // es wird also einfach elementweise multipliziert
  // Am Ende wird diese Matrix mit .matrix() wieder als Matrix interpretiert, also eine DG x NQ matrix
  // und mit phi_gauss.transpose(), einer NQ x 1 - Matrix multipliziert. Im Ergebnis ists dann ein Vektor, also eine DG x 1 - Matrix.
  phiup.row(eid) += dt * (parammap.AdvectionCellTermX[eid].array().rowwise() * vx_gauss.array() +
			  parammap.AdvectionCellTermY[eid].array().rowwise() * vy_gauss.array()).matrix() *
    phi_gauss.transpose();

  
}
////////////////////////////////////////////////// BOUNDARY HANDLING

// could be template over edge... but, sign is changing from left/right, bottom/top
template <int DG>
void boundary_lower(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c, const size_t e)
{
    // GP = DGEDGE
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    // block<1, 2>(e, 0) * PSIe<2,2>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((valueonedge<DG,0>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (-vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 0>;
}
template <int DG>
void boundary_upper(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_X, const size_t c, const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_X.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((valueonedge<DG,2>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 2>;
}
template <int DG>
void boundary_left(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c, const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((valueonedge<DG,3>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (-vel_gauss.array()).max(0));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 3>;
}
template <int DG>
void boundary_right(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup,
    const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c, const size_t e)
{
    LocalEdgeVector<EDGEDOFS(DG)> vel_gauss = normalvel_Y.row(e) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>;
    LocalEdgeVector<EDGEDOFS(DG)> tmp = ((valueonedge<DG,3>(phi, c) * PSIe<EDGEDOFS(DG), EDGEDOFS(DG)>).array() * (vel_gauss.array().max(0)));
    phiup.row(c) -= dt * tmp * PSIe_w<DG, EDGEDOFS(DG), 1>;
}



////////////////////////////////////////////////// EDGE TERMS

/*!
 *
 *  +-----+-----+
 *  |     |     |
 *  |     |     |
 *  |     |     |
 *  +-----+-----+
 *
 *  eid1, eid2 the two elements at the edge
 *  EOE1, EOE2 the position of the edge seen from the element
 *             (0=bottom, 1=right, 2=top, 3=left)
 *  edgeid index of edge
 */
template<int DG>
template<int EOE1, int EOE2>
inline void DGTransport<DG>::edge_term(const ParametricMesh& smesh,
				       const double dt,
				       DGVector<DG>& phiup,
				       const DGVector<DG>& phi, // DG1 (3)
				       const EdgeVector<EDGEDOFS(DG)>& normalvel,
				       const size_t eid1, const size_t eid2,
				       const size_t edgeid)
{
  if (smesh.landmask[eid1]==0) return; // nothing to do if land element
  if (smesh.landmask[eid2]==0) return;
  
  const LocalEdgeVector<GAUSSPOINTS1D(DG)> vel_gauss = normalvel.row(edgeid) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>;

  const LocalEdgeVector<GAUSSPOINTS1D(DG)> tmp =
    (vel_gauss.array().max(0) * (valueonedge<DG,EOE1>(phi, eid1) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>).array() + 
     vel_gauss.array().min(0) * (valueonedge<DG,EOE2>(phi, eid2) * PSIe<EDGEDOFS(DG), GAUSSPOINTS1D(DG)>).array() );
    
  phiup.row(eid1) -= dt * tmp * PSIe_w<DG, GAUSSPOINTS1D(DG), EOE1>;
  phiup.row(eid2) += dt * tmp * PSIe_w<DG, GAUSSPOINTS1D(DG), EOE2>;
}



template <int DG>
void DGTransport<DG>::DGTransportOperator(const ParametricMesh& smesh, const double dt,
    const DGVector<DG>& vx,
    const DGVector<DG>& vy,
    const EdgeVector<EDGEDOFS(DG)>& normalvel_X,
    const EdgeVector<EDGEDOFS(DG)>& normalvel_Y,
    const DGVector<DG>& phi, DGVector<DG>& phiup)
{
    phiup.zero();

    // Cell terms
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        cell_term(smesh, dt, phiup, phi, vx, vy, eid);

    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        size_t ic = iy * smesh.nx; // first index of left cell in row
        size_t ie = iy * (smesh.nx + 1) + 1; // first index of inner velocity in row
	for (size_t i = 0; i < smesh.nx - 1; ++i, ++ic, ++ie)
	  edge_term<1,3>(smesh, dt, phiup, phi, normalvel_Y, ic, ic + 1, ie);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < smesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        size_t ie = ix + smesh.nx; // first index of inner velocity in column
        for (size_t i = 0; i < smesh.ny - 1; ++i, ic += smesh.nx, ie += smesh.nx)
	  edge_term<2,0>(smesh, dt, phiup, phi, normalvel_X, ic, ic + smesh.nx, ie);
    }
   

    // Periodic boundaries: These are like usual edges within the domain
    //#pragma omp parallel for      <<== not sure about parallelization
    for (size_t pc = 0 ; pc<smesh.periodic.size(); ++pc)
      {
	// hm... possible to use eoe1/2 as template parameters???
	if      ( (smesh.periodic[pc].eoe1==0)&&(smesh.periodic[pc].eoe2==0) )
	  edge_term<0,0>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==0)&&(smesh.periodic[pc].eoe2==1) )
	  edge_term<0,1>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==0)&&(smesh.periodic[pc].eoe2==2) )
	  edge_term<0,2>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==0)&&(smesh.periodic[pc].eoe2==3) )
	  edge_term<0,3>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);

	else if ( (smesh.periodic[pc].eoe1==2)&&(smesh.periodic[pc].eoe2==0) )
	  edge_term<2,0>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==2)&&(smesh.periodic[pc].eoe2==1) )
	  edge_term<2,1>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==2)&&(smesh.periodic[pc].eoe2==2) )
	  edge_term<2,2>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==2)&&(smesh.periodic[pc].eoe2==3) )
	  edge_term<2,3>(smesh, dt, phiup, phi, normalvel_X, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);

	else if ( (smesh.periodic[pc].eoe1==1)&&(smesh.periodic[pc].eoe2==0) )
	  edge_term<1,0>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==1)&&(smesh.periodic[pc].eoe2==1) )
	  edge_term<1,1>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==1)&&(smesh.periodic[pc].eoe2==2) )
	  edge_term<1,2>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==1)&&(smesh.periodic[pc].eoe2==3) )
	  edge_term<1,3>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);

	else if ( (smesh.periodic[pc].eoe1==3)&&(smesh.periodic[pc].eoe2==0) )
	  edge_term<3,0>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==3)&&(smesh.periodic[pc].eoe2==1) )
	  edge_term<3,1>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==3)&&(smesh.periodic[pc].eoe2==2) )
	  edge_term<3,2>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else if ( (smesh.periodic[pc].eoe1==3)&&(smesh.periodic[pc].eoe2==3) )
	  edge_term<3,3>(smesh, dt, phiup, phi, normalvel_Y, smesh.periodic[pc].eid1, smesh.periodic[pc].eid2, smesh.periodic[pc].edgeid);
	else
	  abort();
      }

    
    // open boundaries: integrate one-sided edge term
    for (size_t seg = 0 ; seg<4; ++seg) // run over the 4 segments (bot, right, top, left)
      {
#pragma omp parallel for
	for (size_t i = 0; i<smesh.open[seg].size();++i)
	  {
	    const size_t eid = smesh.open[seg][i];
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
		std::cerr << "Wrong open boundary information in the mesh. Boundary side " << smesh.open[seg][i] << " not valid" << std::endl;
		abort();
	      }
	  }
      }


    // Dirichlet: nothing to be done in advection as v=0 on the boundary
   
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
      phiup.row(eid) =  parammap.InverseDGMassMatrix[eid] * phiup.row(eid).transpose();
}

template <int DG>
void DGTransport<DG>::step_rk1(const double dt, DGVector<DG>& phi)
{
    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1);

    phi += tmp1;
}

template <int DG>
void DGTransport<DG>::step_rk2(const double dt, DGVector<DG>& phi)
{
    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp1); // tmp1 = k * F(u)

    phi += tmp1; // phi = phi + k * F(u)     (i.e.: implicit Euler)

    DGTransportOperator(smesh, dt, velx, vely, normalvel_X, normalvel_Y, phi, tmp2); // tmp1 = k * F( u + k * F(u) )

    phi += 0.5 * (tmp2 - tmp1);
}

template <int DG>
void DGTransport<DG>::step_rk3(const double dt, DGVector<DG>& phi)
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

template <int DG>
void DGTransport<DG>::step(const double dt, DGVector<DG>& phi)
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

#undef EDGEDOFS

template class DGTransport<1>;
template class DGTransport<3>;
template class DGTransport<6>;
template class DGTransport<8>;

template void DGTransport<1>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<1>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<3>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<3>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<6>::prepareAdvection(const CGVector<1>& cg_vx, const CGVector<1>& cg_vy);
template void DGTransport<6>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);
template void DGTransport<8>::prepareAdvection(const CGVector<2>& cg_vx, const CGVector<2>& cg_vy);

} /* namespace Nextsim */
