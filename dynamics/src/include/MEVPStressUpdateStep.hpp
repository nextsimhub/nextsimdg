/*!
 * @file MEVPStressUpdateStep.hpp
 *
 * @date Feb 1, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MEVPSTRESSUPDATESTEP_HPP
#define MEVPSTRESSUPDATESTEP_HPP

#include "include/StressUpdateStep.hpp"

#include "include/ParametricMap.hpp"
#include "include/VPParameters.hpp"
#include "include/codeGenerationDGinGauss.hpp"

namespace Nextsim {

template <int DGadvection, int DGstress, int CG>
class MEVPStressUpdateStep : public StressUpdateStep<DGadvection, DGstress> {

public:
    typedef std::array<std::reference_wrapper<DGVector<DGstress>>, N_TENSOR_ELEMENTS>
        SymmetricTensorVector;
    MEVPStressUpdateStep()
        : pmap(nullptr)
    {
    }
    ~MEVPStressUpdateStep() = default;
    void stressUpdateHighOrder(const DynamicsParameters& params, const ParametricMesh& smesh,
        SymmetricTensorVector& stress, const SymmetricTensorVector& strain,
        const DGVector<DGadvection>& H, const DGVector<DGadvection>& A,
        const double deltaT) override
    {
        // Unwrap references
        DGVector<DGstress>& S11 = stress[I11];
        DGVector<DGstress>& S12 = stress[I12];
        DGVector<DGstress>& S22 = stress[I22];

        DGVector<DGstress>& E11 = strain[I11];
        DGVector<DGstress>& E12 = strain[I12];
        DGVector<DGstress>& E22 = strain[I22];

        const VPParameters& vpparameters = reinterpret_cast<const VPParameters&>(params);
        // Number of Gauss points
       constexpr int NGP = ((DGstress == 8) || (DGstress == 6)) ? 3 : (DGstress == 3 ? 2 : -1);
        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            // double DELTA = sqrt(SQR(vpparameters.DeltaMin) + 1.25 * (SQR(E11(i, 0)) + SQR(E22(i,
            // 0)))
            //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
            //   assert(DELTA > 0);

            //   //! Ice strength
            //   double P = vpparameters.Pstar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));

            auto h_gauss = (H.row(i) * PSI<DGadvection, NGP>).array().max(0.0).matrix();
            auto a_gauss = (A.row(i) * PSI<DGadvection, NGP>).array().max(0.0).min(1.0).matrix();

            LocalEdgeVector<NGP* NGP> P
                = (vpparameters.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp())
                      .matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

            // std::cout << P << std::endl << DELTA << std::endl;
            // abort();

            //   // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)

            //  zeta = P / 2.0 / DELTA;
            //  eta = zeta / 4;
            //  (zeta-eta) = zeta - 1/4 zeta = 3/4 zeta

            //   S11.row(i)
            //       += 1.0 / alpha * (2. * eta * E11.row(i) + (zeta - eta) * (E11.row(i) +
            //       E22.row(i)));
            //   S11(i, 0) -= 1.0 / alpha * 0.5 * P;

            // const Eigen::Matrix<Nextsim::FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
            // // get the inverse of the mass matrix scaled with the test-functions in the gauss
            // points,
            // // with the gauss weights and with J. This is a 8 x 9 matrix
            // const Eigen::Matrix<Nextsim::FloatType, 8, 9> imass_psi =
            // ParametricTools::massMatrix<8>(smesh, i).inverse()
            //     * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J.array())).matrix();

            const LocalEdgeVector<NGP* NGP> e11_gauss = E11.row(i) * PSI<DGstress, NGP>;
            const LocalEdgeVector<NGP* NGP> e12_gauss = E12.row(i) * PSI<DGstress, NGP>;
            const LocalEdgeVector<NGP* NGP> e22_gauss = E22.row(i) * PSI<DGstress, NGP>;

            const auto DELTA = (SQR(vpparameters.DeltaMin)
                + 1.25 * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50 * e11_gauss.array() * e22_gauss.array() + e12_gauss.array().square())
                                   .sqrt()
                                   .matrix();

            
            const FloatType alphaInv = 1.0 / alpha;
            const FloatType fac = 1.0 - alphaInv;
            const LocalEdgeVector<NGP* NGP> PDelta = P.array() / DELTA.array();

            S11.row(i) = fac * S11.row(i)
                + (pmap->iMJwPSI[i]
                    * (alphaInv
                        * (PDelta.array()
                                * ((5.0 / 8.0) * e11_gauss.array()
                                    + (3.0 / 8.0) * e22_gauss.array())
                            - 0.5 * P.array())
                              .matrix()
                              .transpose()))
                      .transpose();
            S22.row(i) = fac * S22.row(i)
                + (pmap->iMJwPSI[i]
                    * (alphaInv
                        * (PDelta.array()
                                * ((5.0 / 8.0) * e22_gauss.array()
                                    + (3.0 / 8.0) * e11_gauss.array())
                            - 0.5 * P.array())
                              .matrix()
                              .transpose()))
                      .transpose();
            S12.row(i) = fac * S12.row(i)
                + (pmap->iMJwPSI[i]
                    * (alphaInv
                        * (PDelta.array() * (1.0 / 4.0) * e12_gauss.array()).matrix().transpose()))
                      .transpose();
        }

    }
    void setPMap(ParametricMomentumMap<CG>* pmapIn) { pmap = pmapIn; }

protected:
    ParametricMomentumMap<CG>* pmap;

private:
    //! MEVP parameters
    double alpha = 1500.0;
    double beta = 1500.0;
};

} /* namespace Nextsim */

#endif /* MEVPSTRESSUPDATESTEP_HPP */
