/*!
 * @file mevp.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.no>
 */

#ifndef __MEVP_HPP
#define __MEVP_HPP

#include "dgVector.hpp"
#include "codeGenerationDGinGauss.hpp"

namespace Nextsim {

/*!
 * This namespace collects the routines required for the mEVP solver
 */
namespace mEVP {

    inline constexpr double SQR(double x) { return x * x; }

    template <int DGstress, int DGtracer>
    void StressUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, const double Pstar, const double DeltaMin,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            double DELTA = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            assert(DELTA > 0);

            //! Ice strength
            double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

            double zeta = P / 2.0 / DELTA;
            double eta = zeta / 4;

            // replacement pressure
            P = P * DELTA / (DeltaMin + DELTA);

            // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)
            S11.row(i) *= (1.0 - 1.0 / alpha);
            S12.row(i) *= (1.0 - 1.0 / alpha);
            S22.row(i) *= (1.0 - 1.0 / alpha);

            S11.row(i)
                += 1.0 / alpha * (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            S11(i, 0) -= 1.0 / alpha * 0.5 * P;

            S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));

            S22.row(i)
                += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            S22(i, 0) -= 1.0 / alpha * 0.5 * P;
        }
    }


  //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const Mesh& mesh, CellVector<6>& S11, CellVector<6>& S12,
        CellVector<6>& S22, const CellVector<6>& E11, const CellVector<6>& E12,
        const CellVector<6>& E22, const CellVector<3>& H,
        const CellVector<3>& A, const double Pstar, const double DeltaMin,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {


	  // Here, one should check if it is enough to use a 2-point Gauss rule.
	  // We're dealing with dG2, 3-point Gauss should be required.
	  
	  const LocalEdgeVector<9> h_gauss = (H.block<1,3>(i,0) * BiG33).array().max(0.0).matrix();
	  const LocalEdgeVector<9> a_gauss = (A.block<1,3>(i,0) * BiG33).array().max(0.0).min(1.0).matrix();

	  const LocalEdgeVector<9> e11_gauss = E11.block<1,6>(i,0) * BiG63;
	  const LocalEdgeVector<9> e12_gauss = E12.block<1,6>(i,0) * BiG63;
	  const LocalEdgeVector<9> e22_gauss = E22.block<1,6>(i,0) * BiG63;
	  

	  const LocalEdgeVector<9> DELTA = (SQR(DeltaMin)
	   				    + 1.25* (e11_gauss.array().square() + e22_gauss.array().square())
	   				    + 1.50 * e11_gauss.array() * e22_gauss.array()
	   				    + e12_gauss.array().square()
	   				    ).sqrt().matrix();
	  // double DELTA = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
	  //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
	  //   assert(DELTA > 0);
	  
	  //   //! Ice strength
	  //   double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
	  const LocalEdgeVector<9> P = (Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

          // //   double zeta = P / 2.0 / DELTA;
          // //   double eta = zeta / 4;

          // //   // replacement pressure
          // //   P = P * DELTA / (DeltaMin + DELTA);

	  // std::cout << P << std::endl << DELTA << std::endl;
	  // abort();

          //   // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)
	  S11.row(i) *= (1.0 - 1.0 / alpha);
	  S12.row(i) *= (1.0 - 1.0 / alpha);
	  S22.row(i) *= (1.0 - 1.0 / alpha);

	  //  zeta = P / 2.0 / DELTA;
	  //  eta = zeta / 4;
	  //  (zeta-eta) = zeta - 1/4 zeta = 3/4 zeta
	  
	  //   S11.row(i)
          //       += 1.0 / alpha * (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
          //   S11(i, 0) -= 1.0 / alpha * 0.5 * P;
	  S11.row(i) += 1.0/alpha * ( P.array()/8.0/DELTA.array() * (5.0*e11_gauss.array()+3.0*e22_gauss.array()) - 0.5 * P.array() ).matrix() * IBC63;
	  
          //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
	  // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
	  S12.row(i) += 1.0/alpha * ( P.array()/4.0/DELTA.array() * e12_gauss.array() ).matrix() * IBC63;
	
	  //   S22.row(i)
	  //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
	  //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
	  S22.row(i) += 1.0/alpha * ( P.array()/8.0/DELTA.array() * (5.0*e22_gauss.array()+3.0*e11_gauss.array()) - 0.5 * P.array() ).matrix() * IBC63;
        }
    }



    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const Mesh& mesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<3>& H,
        const CellVector<3>& A, const double Pstar, const double DeltaMin,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {


	  // Here, one should check if it is enough to use a 2-point Gauss rule.
	  // We're dealing with dG2, 3-point Gauss should be required.
	  
	  const LocalEdgeVector<9> h_gauss = (H.block<1,3>(i,0) * BiG33).array().max(0.0).matrix();
	  const LocalEdgeVector<9> a_gauss = (A.block<1,3>(i,0) * BiG33).array().max(0.0).min(1.0).matrix();

	  const LocalEdgeVector<9> e11_gauss = E11.block<1,8>(i,0) * BiG83;
	  const LocalEdgeVector<9> e12_gauss = E12.block<1,8>(i,0) * BiG83;
	  const LocalEdgeVector<9> e22_gauss = E22.block<1,8>(i,0) * BiG83;
	  

	  const LocalEdgeVector<9> DELTA = (SQR(DeltaMin)
	   				    + 1.25* (e11_gauss.array().square() + e22_gauss.array().square())
	   				    + 1.50 * e11_gauss.array() * e22_gauss.array()
	   				    + e12_gauss.array().square()
	   				    ).sqrt().matrix();
	  // double DELTA = sqrt(SQR(DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
	  //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
	  //   assert(DELTA > 0);
	  
	  //   //! Ice strength
	  //   double P = Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
	  const LocalEdgeVector<9> P = (Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

          // //   double zeta = P / 2.0 / DELTA;
          // //   double eta = zeta / 4;

          // //   // replacement pressure
          // //   P = P * DELTA / (DeltaMin + DELTA);

	  // std::cout << P << std::endl << DELTA << std::endl;
	  // abort();

          //   // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)
	  S11.row(i) *= (1.0 - 1.0 / alpha);
	  S12.row(i) *= (1.0 - 1.0 / alpha);
	  S22.row(i) *= (1.0 - 1.0 / alpha);

	  //  zeta = P / 2.0 / DELTA;
	  //  eta = zeta / 4;
	  //  (zeta-eta) = zeta - 1/4 zeta = 3/4 zeta
	  
	  //   S11.row(i)
          //       += 1.0 / alpha * (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
          //   S11(i, 0) -= 1.0 / alpha * 0.5 * P;
	  S11.row(i) += 1.0/alpha * ( P.array()/8.0/DELTA.array() * (5.0*e11_gauss.array()+3.0*e22_gauss.array()) - 0.5 * P.array() ).matrix() * IBC83;
	  
          //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
	  // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
	  S12.row(i) += 1.0/alpha * ( P.array()/4.0/DELTA.array() * e12_gauss.array() ).matrix() * IBC83;
	
	  //   S22.row(i)
	  //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
	  //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
	  S22.row(i) += 1.0/alpha * ( P.array()/8.0/DELTA.array() * (5.0*e22_gauss.array()+3.0*e11_gauss.array()) - 0.5 * P.array() ).matrix() * IBC83;
        }
    }

} /* namespace mEVP */

} /* namespace Nextsim */

#endif /* __MEVP_HPP */
