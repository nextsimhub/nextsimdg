/*----------------------------   mevp.hpp     ---------------------------*/
#ifndef __mevp_HPP
#define __mevp_HPP
/*----------------------------   mevp.hpp     ---------------------------*/

#include "dgvector.hpp"

namespace Nextsim {

  /*!
   * This namespace collects the routines required for the mEVP solver
   */
  namespace mEVP
  {

    inline constexpr double SQR(double x)
    {
      return x * x;
    }
    
    
    template<int DGstress, int DGtracer>
    void StressUpdate(const Mesh& mesh,
		      CellVector<DGstress>& S11, CellVector<DGstress>& S12, CellVector<DGstress>& S22,
		      const CellVector<DGstress>& E11, const CellVector<DGstress>& E12, const CellVector<DGstress>& E22,
		      const CellVector<DGtracer>& H, const CellVector<DGtracer>& A,
		      const double Pstar, const double DeltaMin,
		      const double alpha, const double beta)
    {
      	    
      //! Stress Update
#pragma omp parallel for
      for (size_t i = 0; i < mesh.n; ++i) {
	double DELTA =sqrt(
			   SQR(DeltaMin)
			   + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
			   + 1.50 * E11(i, 0) * E22(i, 0)
			   + SQR(E12(i, 0)));
	assert(DELTA > 0);

	//! Ice strength
	double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
	
	double zeta = P / 2.0 / DELTA;
	double eta = zeta / 4;

	// replacement pressure
	//                P = P * DELTA(i, 0) / (ReferenceScale::DeltaMin + DELTA(i, 0));
	
	//SHEAR(i, 0) = sqrt((SQR(ReferenceScale::DeltaMin) + SQR(E11(i, 0) - E22(i, 0)) + 4.0 * SQR(E12(i, 0))));
	
	// S = S_old + 1/alpha (S(u)-S_old)
	//   = (1-1/alpha) S_old + 1/alpha S(u)
	S11.row(i) *= (1.0 - 1.0 / alpha);
	S12.row(i) *= (1.0 - 1.0 / alpha);
	S22.row(i) *= (1.0 - 1.0 / alpha);
	
	S11.row(i) += 1.0 / alpha * (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
	S11(i, 0) -= 1.0 / alpha * 0.5 * P;
	
	S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
	
	S22.row(i) += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
	S22(i, 0) -= 1.0 / alpha * 0.5 * P;
	
	// // prepare for ellipse-output
	// if ((timestep % NT_vtk == 0)) {
	//   S1(i) = (S11(i, 0) + S22(i, 0)) / P;
	//   S2(i) = sqrt(SQR(S12(i, 0)) + SQR(0.5 * (S11(i, 0) - S22(i, 0)))) / P;
	// }
      }
      
    }
    
  }


}

/*----------------------------   cgmomentum.hpp     ---------------------------*/
/* end of #ifndef __cgmomentum_HPP */
#endif
/*----------------------------   cgmomentum.hpp     ---------------------------*/
