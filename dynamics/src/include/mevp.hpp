/*!
 * @file mevp.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __MEVP_HPP
#define __MEVP_HPP

#include "VPParameters.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

/*!
 * This namespace collects the routines required for the mEVP solver
 */
namespace mEVP {

    inline constexpr double SQR(double x) { return x * x; }

    template <int DGstress, int DGtracer>
    void StressUpdate(const Mesh& mesh,
        const VPParameters& vpparameters,
        CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {
            double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
                + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            assert(DELTA > 0);

            //! Ice strength
            double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

            double zeta = P / 2.0 / DELTA;
            double eta = zeta / 4;

            // replacement pressure
            P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<6>& S11, CellVector<6>& S12,
        CellVector<6>& S22, const CellVector<6>& E11, const CellVector<6>& E12,
        const CellVector<6>& E22, const CellVector<3>& H,
        const CellVector<3>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 6>(i, 0) * PSI<6,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 6>(i, 0) * PSI<6,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 6>(i, 0) * PSI<6,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC63 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC63 * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC63 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix()).transpose();
        }
    }

    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<3>& S11, CellVector<3>& S12,
        CellVector<3>& S22, const CellVector<3>& E11, const CellVector<3>& E12,
        const CellVector<3>& E22, const CellVector<3>& H,
        const CellVector<3>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 3>(i, 0) * PSI<3,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 3>(i, 0) * PSI<3,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 3>(i, 0) * PSI<3,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();

            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC33 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC33 * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC33 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix()).transpose();
        }
    }

    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<3>& S11, CellVector<3>& S12,
        CellVector<3>& S22, const CellVector<3>& E11, const CellVector<3>& E12,
        const CellVector<3>& E22, const CellVector<1>& H,
        const CellVector<1>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG???, less points???

            const LocalEdgeVector<9> h_gauss = (H(i, 0) * PSI<1,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A(i, 0) * PSI<1,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 3>(i, 0) * PSI<3,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 3>(i, 0) * PSI<3,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 3>(i, 0) * PSI<3,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();

            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC33 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC33 * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC33 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix()).transpose();
        }
    }

    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<1>& S11, CellVector<1>& S12,
        CellVector<1>& S22, const CellVector<1>& E11, const CellVector<1>& E12,
        const CellVector<1>& E22, const CellVector<3>& H,
        const CellVector<3>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            LocalEdgeVector<9> ones;
            ones << 1., 1., 1., 1., 1., 1., 1., 1., 1.;
            const LocalEdgeVector<9> e11_gauss = E11(i, 0) * ones;
            const LocalEdgeVector<9> e12_gauss = E12(i, 0) * ones;
            const LocalEdgeVector<9> e22_gauss = E22(i, 0) * ones;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();

            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC13 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC13 * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC13 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix()).transpose();
        }
    }

    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<3>& H,
        const CellVector<3>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar
                * h_gauss.array()
                * (-20.0 * (1.0 - a_gauss.array())).exp())
                                             .matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC83 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC83 * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC83 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix()).transpose();
        }
    }

    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<6>& H,
        const CellVector<6>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 6>(i, 0) * PSI<6,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 6>(i, 0) * PSI<6,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC83 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC83 * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC83 * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix()).transpose();
        }
    }

    //! The main EVP stress update using higher order represenatation of A and H
    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const Mesh& mesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<1>& H,
        const CellVector<1>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const double h_gauss = H(i, 0);
            const double a_gauss = A(i, 0);

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const double P = vpparameters.Pstar * h_gauss * exp(-20.0 * (1.0 - a_gauss));

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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
            S11.row(i) += IBC83 * (1.0 / alpha * (P / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P).matrix()).transpose();

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += IBC83 * (1.0 / alpha * (P / 4.0 / DELTA.array() * e12_gauss.array()).matrix()).transpose();

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += IBC83 * (1.0 / alpha * (P / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P).matrix()).transpose();
        }
    }

    // Stress Update (SasipMesh)

    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const ParametricTransformation<2, 8>& ptrans,
        const SasipMesh& smesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<3>& H,
        const CellVector<3>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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

            // const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
            // // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // // with the gauss weights and with J. This is a 8 x 9 matrix
            // const Eigen::Matrix<Nextsim::FloatType, 8, 9> imass_psi = ParametricTools::massMatrix<8>(smesh, i).inverse()
            //     * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS_3.array() * J.array())).matrix();

            S11.row(i) += ptrans.iMJwPSI[i] * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix().transpose());

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += ptrans.iMJwPSI[i] * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix().transpose());

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += ptrans.iMJwPSI[i] * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix().transpose());
        }
    }

    void StressUpdateHighOrder(const VPParameters& vpparameters,
        const SasipMesh& smesh, CellVector<8>& S11, CellVector<8>& S12,
        CellVector<8>& S22, const CellVector<8>& E11, const CellVector<8>& E12,
        const CellVector<8>& E22, const CellVector<3>& H,
        const CellVector<3>& A,
        const double alpha, const double beta)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<9> h_gauss = (H.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).matrix();
            const LocalEdgeVector<9> a_gauss = (A.block<1, 3>(i, 0) * PSI<3,3>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<9> e11_gauss = E11.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e12_gauss = E12.block<1, 8>(i, 0) * PSI<8,3>;
            const LocalEdgeVector<9> e22_gauss = E22.block<1, 8>(i, 0) * PSI<8,3>;

            const LocalEdgeVector<9> DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array()
                + e12_gauss.array().square())
                                                 .sqrt()
                                                 .matrix();
            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i, 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
            const LocalEdgeVector<9> P = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

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

            const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
            // get the inverse of the mass matrix scaled with the test-functions in the gauss points,
            // with the gauss weights and with J. This is a 8 x 9 matrix
            const Eigen::Matrix<Nextsim::FloatType, 8, 9> imass_psi = ParametricTools::massMatrix<8>(smesh, i).inverse()
                * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS_3.array() * J.array())).matrix();

            S11.row(i) += imass_psi * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix().transpose());

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            S12.row(i) += imass_psi * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix().transpose());

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            S22.row(i) += imass_psi * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix().transpose());
        }
    }

} /* namespace mEVP */

} /* namespace Nextsim */

#endif /* __MEVP_HPP */
