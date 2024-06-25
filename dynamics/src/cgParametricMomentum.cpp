/*!
 * @file ParametricMomentum.cpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "cgParametricMomentum.hpp"
#include "Interpolations.hpp"
#include "MEB.hpp"
#include "BBM.hpp"
#include "ParametricTools.hpp"
#include "mevp.hpp"
#include "dgVisu.hpp"
#include "VectorManipulations.hpp"

namespace Nextsim {

#define DGSTRESS(CG) ( (CG==1?3:(CG==2?8:-1) ) )

  ////////////////////////////////////////////////// Strain (ParametricMesh)

  template <int CG>
  void CGParametricMomentum<CG>::ProjectCGVelocityToDGStrain()
  {
    // !!! must still be converted to the spherical system!!!
  
    assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == vx.rows());
    assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == vy.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E11.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E12.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E22.rows());

    const int cgshift = CG * smesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < smesh.ny; ++row) { 
      int dgi = smesh.nx * row; //!< Index of dg vector
      int cgi = CG * cgshift * row; //!< Lower left index of cg vector

      for (size_t col = 0; col < smesh.nx; ++col, ++dgi, cgi += CG) { // loop over all elements

	if (smesh.landmask[dgi]==0) // only on ice
	  continue;
	      
	      
	// get the 4 (cg1) 9 (cg2) local x/y - velocity coefficients on the element
	Eigen::Matrix<double, CGDOFS(CG), 1> vx_local, vy_local;
	if (CG == 1) {
	  vx_local << vx(cgi), vx(cgi + 1), vx(cgi + cgshift), vx(cgi + 1 + cgshift);
	  vy_local << vy(cgi), vy(cgi + 1), vy(cgi + cgshift), vy(cgi + 1 + cgshift);
	} else if (CG == 2) {
	  vx_local << vx(cgi), vx(cgi + 1), vx(cgi + 2), vx(cgi + cgshift), vx(cgi + 1 + cgshift),
	    vx(cgi + 2 + cgshift), vx(cgi + 2 * cgshift), vx(cgi + 1 + 2 * cgshift),
	    vx(cgi + 2 + 2 * cgshift);

	  vy_local << vy(cgi), vy(cgi + 1), vy(cgi + 2), vy(cgi + cgshift), vy(cgi + 1 + cgshift),
	    vy(cgi + 2 + cgshift), vy(cgi + 2 * cgshift), vy(cgi + 1 + 2 * cgshift),
	    vy(cgi + 2 + 2 * cgshift);
	} else
	  abort();

	// Solve (E, Psi) = (0.5(DV + DV^T), Psi)
	// by integrating rhs and inverting with dG(stress) mass matrix
	//
	E11.row(dgi) = pmap.iMgradX[dgi] * vx_local;
	E22.row(dgi) = pmap.iMgradY[dgi] * vy_local;
	E12.row(dgi) = 0.5 * (pmap.iMgradX[dgi] * vy_local + pmap.iMgradY[dgi] * vx_local);

	if (smesh.CoordinateSystem == SPHERICAL)
	  {
	    E11.row(dgi) -= pmap.iMM[dgi] * vy_local;
	    E12.row(dgi) += 0.5 * pmap.iMM[dgi] * vx_local;
	  }
      }
    }
  }

  ////////////////////////////////////////////////// STRESS Tensor
  // Sasip-Mesh Interface
  template <int CG>
  void CGParametricMomentum<CG>::DivergenceOfStress(const double scale, CGVector<CG>& tx,
						    CGVector<CG>& ty) const
  {
#pragma omp parallel for
    for (size_t i=0;i<tx.rows();++i)
      {
	tx(i)=0.0;
	ty(i)=0.0;
      }
    
    // parallelization in stripes
    for (size_t p = 0; p < 2; ++p)
#pragma omp parallel for schedule(static)
      for (size_t cy = 0; cy < smesh.ny; ++cy) //!< loop over all cells of the mesh
        {
	  if (cy % 2 == p) {
	    size_t c = smesh.nx * cy;
	    for (size_t cx = 0; cx < smesh.nx; ++cx, ++c) //!< loop over all cells of the mesh
	      if (smesh.landmask[c]==1) // only on ice!
		AddStressTensorCell(scale, c, cx, cy, tx, ty);
	  }
        }
    // set zero on the Dirichlet boundaries
    DirichletZero(tx);
    DirichletZero(ty);
    FreeSlipZero(tx,ty);
    
    // add the contributions on the periodic boundaries
    VectorManipulations::CGAveragePeriodic(smesh, tx);
    VectorManipulations::CGAveragePeriodic(smesh, ty);
  }

  template <int CG>
  void CGParametricMomentum<CG>::DirichletZero(CGVector<CG>& v) const
  {
    // the four segments bottom, right, top, left, are each processed in parallel
    for (size_t seg=0;seg<4;++seg)
      {
#pragma omp parallel for
        for (size_t i = 0; i < smesh.dirichlet[seg].size(); ++i) {

	  const size_t eid = smesh.dirichlet[seg][i];
	  const size_t ix = eid % smesh.nx; // compute 'coordinate' of element
	  const size_t iy = eid / smesh.nx;

	  if (seg == 0) // bottom
	    for (size_t j = 0; j < CG + 1; ++j)
	      v(iy * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) = 0.0;
	  else if (seg == 1) // right
	    for (size_t j = 0; j < CG + 1; ++j)
	      v(iy * CG * (CG * smesh.nx + 1) + CG * ix + CG + (CG * smesh.nx + 1) * j, 0) = 0.0;
	  else if (seg == 2) // top
	    for (size_t j = 0; j < CG + 1; ++j)
	      v((iy + 1) * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) = 0.0;
	  else if (seg == 3) // left
	    for (size_t j = 0; j < CG + 1; ++j)
	      v(iy * CG * (CG * smesh.nx + 1) + CG * ix + (CG * smesh.nx + 1) * j, 0) = 0.0;
	  else {
	    std::cerr << "That should not have happened!" << std::endl;
	    abort();
	  }
        }
      }
  }

  
  template <int CG>
  void CGParametricMomentum<CG>::FreeSlipZero(CGVector<CG>& vx,CGVector<CG>& vy) const
  {
    // the four segments bottom, right, top, left, are each processed in parallel
    for (size_t seg=0;seg<4;++seg)
      {
#pragma omp parallel for
        for (size_t i = 0; i < smesh.freeslip[seg].size(); ++i) {

	  const size_t eid = smesh.freeslip[seg][i];
	  const size_t ix = eid % smesh.nx; // compute 'coordinate' of element
	  const size_t iy = eid / smesh.nx;


	  for (size_t j = 0; j < CG + 1; ++j)
	    {
	      // normal vector on edge in dof (j)
	      const Eigen::Matrix<Nextsim::FloatType, 1, 2> normal = smesh.cgnormal<CG>(eid,seg,j);

	      if (seg == 0) // bottom
		{
		  const double vn =
		    normal(0,0) * vx(iy * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) +
		    normal(0,1) * vy(iy * CG * (CG * smesh.nx + 1) + CG * ix + j, 0);
		  vx(iy * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) -= vn * normal(0,0);
		  vy(iy * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) -= vn * normal(0,1);
		}
	      else if (seg == 1) // right
		{
		  const double vn =
		    normal(0,0) * vx(iy * CG * (CG * smesh.nx + 1) + CG * ix + CG + (CG * smesh.nx + 1) * j, 0) +
		    normal(0,1) * vy(iy * CG * (CG * smesh.nx + 1) + CG * ix + CG + (CG * smesh.nx + 1) * j, 0);
		  vx(iy * CG * (CG * smesh.nx + 1) + CG * ix + CG + (CG * smesh.nx + 1) * j, 0) -= vn * normal(0,0);
		  vy(iy * CG * (CG * smesh.nx + 1) + CG * ix + CG + (CG * smesh.nx + 1) * j, 0) -= vn * normal(0,1);
		}
	      else if (seg == 2) // top
		{
		  const double vn =
		    normal(0,0) * vx((iy + 1) * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) +
		    normal(0,1) * vy((iy + 1) * CG * (CG * smesh.nx + 1) + CG * ix + j, 0);
		  vx((iy + 1) * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) -= vn * normal(0,0);
		  vy((iy + 1) * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) -= vn * normal(0,1);
		}
	      else if (seg == 3) // left
		{
		  const double vn = 
		    normal(0,0) * vx(iy * CG * (CG * smesh.nx + 1) + CG * ix + (CG * smesh.nx + 1) * j, 0) +
		    normal(0,1) * vy(iy * CG * (CG * smesh.nx + 1) + CG * ix + (CG * smesh.nx + 1) * j, 0);
		  vx(iy * CG * (CG * smesh.nx + 1) + CG * ix + (CG * smesh.nx + 1) * j, 0) -= vn * normal(0,0);
		  vy(iy * CG * (CG * smesh.nx + 1) + CG * ix + (CG * smesh.nx + 1) * j, 0) -= vn * normal(0,1);
		}
	      else {
		std::cerr << "That should not have happened!" << std::endl;
		abort();
	      }
	    }
	}
      }
  }


// --------------------------------------------------

  template <int CG>
  template <int DG>
  void CGParametricMomentum<CG>::prepareIteration(const DGVector<DG>& H, const DGVector<DG>& A)
  {
    // copy old velocity
    vx_mevp = vx;
    vy_mevp = vy;
    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(smesh, cg_A, A);
    VectorManipulations::CGAveragePeriodic(smesh, cg_A);
    Interpolations::DG2CG(smesh, cg_H, H);
    VectorManipulations::CGAveragePeriodic(smesh, cg_H);

    // limit A to [0,1] and H to [0, ...)
    cg_A = cg_A.cwiseMin(1.0);
    cg_A = cg_A.cwiseMax(1.e-4);
    cg_H = cg_H.cwiseMax(1.e-4);
}

  template <int CG>
  template <int DG>
  void CGParametricMomentum<CG>::mEVPStep(const VPParameters& params,
					  const size_t NT_evp, const double alpha, const double beta,
					  double dt_adv,
					  const DGVector<DG>& H, const DGVector<DG>& A)
  {
    
    // Compute Strain Rate
    
    ProjectCGVelocityToDGStrain();
    


    // Update the stresses according to the mEVP model
    
    Nextsim::mEVP::StressUpdateHighOrder(params, pmap, smesh, S11, S12, S22, E11, E12, E22, H, A, alpha, beta);
    

    // Compute the divergence of the stress tensor
    

    double stressscale = 1.0; // 2nd-order Stress term has different scaling with the EarthRadius
    //    if (smesh.CoordinateSystem == Nextsim::SPHERICAL)
    //      stressscale = 1.0/Nextsim::EarthRadius/Nextsim::EarthRadius;

    DivergenceOfStress(stressscale, tmpx, tmpy); // Compute divergence of stress tensor
    
    

    // Update the velocity

    //	    update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
    for (int i = 0; i < vx.rows(); ++i) {
      double absatm = sqrt(ax(i)*ax(i)+ay(i)*ay(i));
      double absocn = sqrt(SQR(vx(i)-ox(i)) + SQR(vy(i)-oy(i)));

      vx(i) = (1.0
	       / (params.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
		  + cg_A(i) * params.F_ocean 
		  * absocn ) // implicit parts
	       * (params.rho_ice * cg_H(i) / dt_adv * (beta * vx(i) + vx_mevp(i)) // pseudo-timestepping
              + cg_A(i)
                  * (params.F_atm * absatm * ax(i) + // atm forcing
		     params.F_ocean * absocn 
		     * ox(i)) // ocean forcing
		  + params.rho_ice * cg_H(i) * params.fc
		  * (vy(i) - oy(i)) // cor + surface
		  + tmpx(i)/pmap.lumpedcgmass(i)
		  ));
      vy(i) = (1.0
	       / (params.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
		  + cg_A(i) * params.F_ocean 
		  * absocn ) // implicit parts
	       * (params.rho_ice * cg_H(i) / dt_adv * (beta * vy(i) + vy_mevp(i)) // pseudo-timestepping
              + cg_A(i)
                  * (params.F_atm * absatm * ay(i) + // atm forcing
		     params.F_ocean * absocn 
		     * oy(i)) // ocean forcing
		  + params.rho_ice * cg_H(i) * params.fc
		  * (ox(i) - vx(i))
		  + tmpy(i)/pmap.lumpedcgmass(i)
		  ));
    }
       
       
    DirichletZero(vx);
    DirichletZero(vy);
    FreeSlipZero(vx,vy);
    
    // Landmask
    const size_t inrow = CG*smesh.nx+1;
#pragma omp parallel for
    for (size_t eid=0;eid<smesh.nelements;++eid)
      if (smesh.landmask[eid]==0)
	{
	  const size_t ex = eid%smesh.nx;
	  const size_t ey = eid/smesh.nx;
	  for (int jy=0;jy<CG+1;++jy)
	    for (int jx=0;jx<CG+1;++jx)
	      {
		vx(inrow*(CG*ey+jy)+CG*ex+jx)=0.0;
		vy(inrow*(CG*ey+jy)+CG*ex+jx)=0.0;
	      }
	}
    
}
  // --------------------------------------------------
  template <int CG>
  template <int DG>
  void CGParametricMomentum<CG>::prepareIteration(const DGVector<DG>& H,
						  const DGVector<DG>& A, const DGVector<DG>& D)
  {

    // set the average sub-iteration velocity to zero
    avg_vx.setZero();
    avg_vy.setZero();   

    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(smesh, cg_A, A);
    VectorManipulations::CGAveragePeriodic(smesh,cg_A);
    Interpolations::DG2CG(smesh, cg_H, H);
    VectorManipulations::CGAveragePeriodic(smesh, cg_H);
    Interpolations::DG2CG(smesh, cg_D, D);
    VectorManipulations::CGAveragePeriodic(smesh, cg_D);

    // limit A and D to [0,1] and H to [0, ...)
    cg_A = cg_A.cwiseMin(1.0);
    cg_A = cg_A.cwiseMax(1.e-4);
    cg_H = cg_H.cwiseMax(1.e-4);
    cg_D = cg_D.cwiseMin(1.0);
    cg_D = cg_D.cwiseMax(0.0);
  }

/* This is Hunke and Dukowicz's solution to (22), multiplied
 * with (\Delta t/m)^2 to ensure stability for c' = 0 
 *
 * This scheme includes Coriolis terms in an implicit way
 */
  template <int CG>
  template <int DG>
  void CGParametricMomentum<CG>::BBMStep(const MEBParameters& params, const size_t NT_meb,
      double dt_adv, const DGVector<DG>& H, const DGVector<DG>& A, DGVector<DG>& D)
  {

    double dt_mom = dt_adv / NT_meb;

    
    //! Compute Strain Rate
    ProjectCGVelocityToDGStrain();
    

    
    // TODO compute stress update with precomputed transformations
    Nextsim::MEB::StressUpdateHighOrder<CG, DGSTRESS(CG), DG>(params, smesh, S11, S12, S22, E11, E12, E22, H, A, D, dt_mom);
    // Nextsim::MEB::StressUpdateHighOrder(params, ptrans, smesh, S11, S12, S22, E11, E12, E22, H, A, D, dt_mom);
    

    
    // Divergence of the stress tensor
#pragma omp parallel for
    for (int i = 0; i < tmpx.rows(); ++i)
        tmpx(i) = tmpy(i) = 0;

    // AddStressTensor(ptrans_stress, -1.0, tmpx, tmpy, S11, S12, S22);
    // AddStressTensor(-1.0, tmpx, tmpy);
    // Check first parameter, this scaling 
    DivergenceOfStress(1.0, tmpx, tmpy); // Compute divergence of stress tensor

    // FIXME: We're missing the gradient of the sea-surface slope (\nabla \eta)


    /* This is Hunke and Dukowicz's solution to (22), multiplied
     * with (\Delta t/m)^2 to ensure stability for c' = 0 */
    double const cos_ocean_turning_angle = std::cos(params.ocean_turning_angle * M_PI / 180.);
    double const sin_ocean_turning_angle = std::sin(params.ocean_turning_angle * M_PI / 180.);

#pragma omp parallel for
    for (int i = 0; i < vx.rows(); ++i) {
        // FIXME: dte_over_mass should include snow (total mass)
        double const dte_over_mass = dt_mom / (params.rho_ice * cg_H(i));
        double const uice = vx(i);
        double const vice = vy(i);

        double const c_prime = cg_A(i) * params.F_ocean * std::hypot(ox(i) - uice, oy(i) - vice);

        // FIXME: Need the grounding term: tau_b = C_bu[i]/(std::hypot(uice,vice)+u0);
        double const tau_b = 0.;
        double const alpha = 1. + dte_over_mass * (c_prime * cos_ocean_turning_angle + tau_b);
        /* FIXME: We need latitude here. Then this becomes:
         * double const beta   = dt_mom*params.fc +
         * dte_over_mass*c_prime*std::copysign(sin_ocean_turning_angle, lat[i]); */
        double const beta = dt_mom * params.fc + dte_over_mass * c_prime * sin_ocean_turning_angle;
        double const rdenom = 1. / (alpha * alpha + beta * beta);

        double const drag_atm = cg_A(i) * params.F_atm * std::hypot(ax(i), ay(i));
        double const tau_x = drag_atm * ax(i)
            + c_prime * (ox(i) * cos_ocean_turning_angle - oy(i) * sin_ocean_turning_angle);
        /* FIXME: Need latitude here. Then This becomes:
         * + c_prime*( ox(i)*cos_ocean_turning_angle - oy(i)*std::copysign(sin_ocean_turning_angle,
         * lat[i]) ); */
        double const tau_y = drag_atm * ay(i)
            + c_prime * (oy(i) * cos_ocean_turning_angle + ox(i) * sin_ocean_turning_angle);
        /* FIXME: Need latitude here. Then This becomes:
         * + c_prime*( oy(i)*cos_ocean_turning_angle + ox(i)*std::copysign(sin_ocean_turning_angle,
         * lat[i]) ); */

        // We need to divide the gradient terms with the lumped mass matrix term
        double const grad_x = tmpx(i) / pmap.lumpedcgmass(i);
        double const grad_y = tmpy(i) / pmap.lumpedcgmass(i);

        vx(i) = alpha * uice + beta * vice
            + dte_over_mass * (alpha * (grad_x + tau_x) + beta * (grad_y + tau_y));
        vx(i) *= rdenom;

        vy(i) = alpha * vice - beta * uice
            + dte_over_mass * (alpha * (grad_y + tau_y) + beta * (grad_x + tau_x));
        vy(i) *= rdenom;
    }

    DirichletZero(vx);
    DirichletZero(vy);
    FreeSlipZero(vx,vy);
    

    avg_vx += vx/NT_meb;
    avg_vy += vy/NT_meb;
}


template <int CG>
template <int DG>
void CGParametricMomentum<CG>::MEBStep(const MEBParameters& params,
    const size_t NT_meb, double dt_adv, const DGVector<DG>& H, const DGVector<DG>& A,
    DGVector<DG>& D)
{

    double dt_mom = dt_adv / NT_meb;

    
    //! Compute Strain Rate
    ProjectCGVelocityToDGStrain();
    

    
    // TODO compute stress update with precomputed transformations
    Nextsim::MEB::StressUpdateHighOrder<CG, DGSTRESS(CG), DG>(params, smesh, S11, S12, S22, E11, E12, E22, H, A, D, dt_mom);
    //Nextsim::BBM::StressUpdateHighOrder<CG, DGSTRESS(CG), DG>(params, smesh, S11, S12, S22, E11, E12, E22, H, A, D, dt_mom);

    // Nextsim::MEB::StressUpdateHighOrder(params, ptrans, smesh, S11, S12, S22, E11, E12, E22, H, A, D, dt_mom);
    

#pragma omp parallel for
    for (int i = 0; i < tmpx.rows(); ++i)
        tmpx(i) = tmpy(i) = 0;

    
    //AddStressTensor(-1.0, tmpx, tmpy);
    DivergenceOfStress(1.0, tmpx, tmpy); // Compute divergence of stress tensor
    


    double const cos_ocean_turning_angle = std::cos(params.ocean_turning_angle * M_PI / 180.);
    double const sin_ocean_turning_angle = std::sin(params.ocean_turning_angle * M_PI / 180.);


    
#pragma omp parallel for
    for (int i = 0; i < vx.rows(); ++i) {

        double absatm = sqrt(ax(i)*ax(i)+ay(i)*ay(i));

        double corsurf_x = vx(i) - ox(i);
        double corsurf_y = vy(i) - oy(i);
        double absocn = sqrt(SQR( corsurf_x ) + SQR( corsurf_y ));

        vx(i) = (1.0
            / (params.rho_ice * cg_H(i) / dt_mom // implicit parts
                + cg_A(i) * params.F_ocean
                    * absocn ) // implicit parts
            * (params.rho_ice * cg_H(i) / dt_mom * vx(i)
                + cg_A(i) * (params.F_atm * absatm * ax(i) + // atm forcing
                      params.F_ocean * absocn * ox(i)) // ocean forcing
                + params.rho_ice * cg_H(i) * params.fc
                    * corsurf_y  )); // cor + surface

        vy(i) = (1.0
            / (params.rho_ice * cg_H(i) / dt_mom // implicit parts
                + cg_A(i) * params.F_ocean
                    * absocn ) // implicit parts
            * (params.rho_ice * cg_H(i) / dt_mom * vy(i)
                + cg_A(i) * (params.F_atm * absatm * ay(i) + // atm forcing
                      params.F_ocean * absocn * oy(i)) // ocean forcing
                + params.rho_ice * cg_H(i) * params.fc
                    * corsurf_x )); // cor + surface

        vx(i) += (1.0
                     / (params.rho_ice * cg_H(i) / dt_mom // implicit parts
                         + cg_A(i) * params.F_ocean
                             * absocn ) // implicit parts
                     * tmpx(i))
            / pmap.lumpedcgmass(i);
        ;

        vy(i) += (1.0
                     / (params.rho_ice * cg_H(i) / dt_mom // implicit parts
                         + cg_A(i) * params.F_ocean
                             * absocn ) // implicit parts
                     * tmpy(i))
            / pmap.lumpedcgmass(i);
        ;
    }
    

    
    DirichletZero(vx);
    DirichletZero(vy);
    FreeSlipZero(vx,vy);

    avg_vx += vx/NT_meb;
    avg_vy += vy/NT_meb;   
  }
  // --------------------------------------------------

  template class CGParametricMomentum<1>;
  template class CGParametricMomentum<2>;

  template void CGParametricMomentum<1>::prepareIteration(const DGVector<1>& H, const DGVector<1>& A);
  template void CGParametricMomentum<1>::prepareIteration(const DGVector<3>& H, const DGVector<3>& A);
  template void CGParametricMomentum<1>::prepareIteration(const DGVector<6>& H, const DGVector<6>& A);
  template void CGParametricMomentum<2>::prepareIteration(const DGVector<1>& H, const DGVector<1>& A);
  template void CGParametricMomentum<2>::prepareIteration(const DGVector<3>& H, const DGVector<3>& A);
  template void CGParametricMomentum<2>::prepareIteration(const DGVector<6>& H, const DGVector<6>& A);

  template void CGParametricMomentum<1>::prepareIteration(const DGVector<1>& H, const DGVector<1>& A, const DGVector<1>& D);
  template void CGParametricMomentum<1>::prepareIteration(const DGVector<3>& H, const DGVector<3>& A, const DGVector<3>& D);
  template void CGParametricMomentum<1>::prepareIteration(const DGVector<6>& H, const DGVector<6>& A, const DGVector<6>& D);
  template void CGParametricMomentum<2>::prepareIteration(const DGVector<1>& H, const DGVector<1>& A, const DGVector<1>& D);
  template void CGParametricMomentum<2>::prepareIteration(const DGVector<3>& H, const DGVector<3>& A, const DGVector<3>& D);
  template void CGParametricMomentum<2>::prepareIteration(const DGVector<6>& H, const DGVector<6>& A, const DGVector<6>& D);

  // --------------------------------------------------

  template void CGParametricMomentum<1>::mEVPStep(const VPParameters& params,
						  size_t NT_evp, double alpha, double beta,
						  double dt_adv,
						  const DGVector<1>& H, const DGVector<1>& A);
  template void CGParametricMomentum<1>::mEVPStep(const VPParameters& params,
						  size_t NT_evp, double alpha, double beta,
						  double dt_adv,
						  const DGVector<3>& H, const DGVector<3>& A);
  template void CGParametricMomentum<1>::mEVPStep(const VPParameters& params,
						  size_t NT_evp, double alpha, double beta,
						  double dt_adv,
						  const DGVector<6>& H, const DGVector<6>& A);

  template void CGParametricMomentum<2>::mEVPStep(const VPParameters& params,
						  size_t NT_evp, double alpha, double beta,
						  double dt_adv,
						  const DGVector<1>& H, const DGVector<1>& A);
  template void CGParametricMomentum<2>::mEVPStep(const VPParameters& params,
						  size_t NT_evp, double alpha, double beta,
						  double dt_adv,
						  const DGVector<3>& H, const DGVector<3>& A);
  template void CGParametricMomentum<2>::mEVPStep(const VPParameters& params,
						  size_t NT_evp, double alpha, double beta,
						  double dt_adv,
						  const DGVector<6>& H, const DGVector<6>& A);

  // --------------------------------------------------

  template void CGParametricMomentum<1>::MEBStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<1>& H, const DGVector<1>& A, DGVector<1>& D);
  template void CGParametricMomentum<1>::MEBStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<3>& H, const DGVector<3>& A, DGVector<3>& D);
  template void CGParametricMomentum<1>::MEBStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<6>& H, const DGVector<6>& A, DGVector<6>& D);


  template void CGParametricMomentum<2>::MEBStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<1>& H, const DGVector<1>& A, DGVector<1>& D);
  template void CGParametricMomentum<2>::MEBStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<3>& H, const DGVector<3>& A, DGVector<3>& D);
  template void CGParametricMomentum<2>::MEBStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<6>& H, const DGVector<6>& A, DGVector<6>& D);

// --------------------------------------------------

  template void CGParametricMomentum<1>::BBMStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<1>& H, const DGVector<1>& A, DGVector<1>& D);
  template void CGParametricMomentum<1>::BBMStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<3>& H, const DGVector<3>& A, DGVector<3>& D);
  template void CGParametricMomentum<1>::BBMStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<6>& H, const DGVector<6>& A, DGVector<6>& D);


  template void CGParametricMomentum<2>::BBMStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<1>& H, const DGVector<1>& A, DGVector<1>& D);
  template void CGParametricMomentum<2>::BBMStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<3>& H, const DGVector<3>& A, DGVector<3>& D);
  template void CGParametricMomentum<2>::BBMStep(const MEBParameters& params,
						 size_t NT_evp, double dt_adv,
						 const DGVector<6>& H, const DGVector<6>& A, DGVector<6>& D);

} /* namespace Nextsim */
