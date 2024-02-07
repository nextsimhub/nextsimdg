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

namespace Nextsim {

template <int DGadvection, int DGstress, int CG>
class MEVPStressUpdateStep : public StressUpdateStep<DGadvection, DGstress> {
//using StressUpdateStep<DGadvection, DGstress>::SymmetricTensorVector;

public:
    typedef std::array<DGVector<DGstress>, N_TENSOR_ELEMENTS> SymmetricTensorVector;
    MEVPStressUpdateStep()
        : pmap (nullptr)
    {
    }
    ~MEVPStressUpdateStep() = default;
    void stressUpdateHighOrder(const DynamicsParameters& params,
            const ParametricMesh& smesh,
            SymmetricTensorVector& stress, const SymmetricTensorVector& strain,
            const DGVector<DGadvection>& h, const DGVector<DGadvection>& a,
            const double deltaT) override
    {
        const VPParameters& vpParams = reinterpret_cast<const VPParameters&>(params);
        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            // Here, one should check if it is enough to use a 2-point Gauss rule.
            // We're dealing with dG2, 3-point Gauss should be required.

            const LocalEdgeVector<nGauss * nGauss> h_gauss = (h.row(i) * PSI<DGadvection, nGauss>).array().max(0.0).matrix();
            const LocalEdgeVector<nGauss * nGauss> a_gauss = (a.row(i) * PSI<DGadvection, nGauss>).array().max(0.0).min(1.0).matrix();

            const LocalEdgeVector<nGauss * nGauss> e11_gauss = strain[I11].row(i) * PSI<DGstress, nGauss>;
            const LocalEdgeVector<nGauss * nGauss> e12_gauss = strain[I12].row(i) * PSI<DGstress, nGauss>;
            const LocalEdgeVector<nGauss * nGauss> e22_gauss = strain[I22].row(i) * PSI<DGstress, nGauss>;

            const LocalEdgeVector<nGauss * nGauss> DELTA = (SQR(vpParams.DeltaMin)
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
            const LocalEdgeVector<nGauss * nGauss> P = (vpParams.Pstar * h_gauss.array() * (-20.0 * (1.0 - a_gauss.array())).exp()).matrix();

            // //   double zeta = P / 2.0 / DELTA;
            // //   double eta = zeta / 4;

            // //   // replacement pressure
            // //   P = P * DELTA / (vpparameters.DeltaMin + DELTA);

            // std::cout << P << std::endl << DELTA << std::endl;
            // abort();

            //   // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)
            stress[I11].row(i) *= (1.0 - 1.0 / alpha);
            stress[I12].row(i) *= (1.0 - 1.0 / alpha);
            stress[I22].row(i) *= (1.0 - 1.0 / alpha);

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
            //     * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J.array())).matrix();

            stress[I11].row(i) += pmap->iMJwPSI[i] * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e11_gauss.array() + 3.0 * e22_gauss.array()) - 0.5 * P.array()).matrix().transpose());

            //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
            // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
            stress[I12].row(i) += pmap->iMJwPSI[i] * (1.0 / alpha * (P.array() / 4.0 / DELTA.array() * e12_gauss.array()).matrix().transpose());

            //   S22.row(i)
            //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) + E22.row(i)));
            //   S22(i, 0) -= 1.0 / alpha * 0.5 * P;
            stress[I22].row(i) += pmap->iMJwPSI[i] * (1.0 / alpha * (P.array() / 8.0 / DELTA.array() * (5.0 * e22_gauss.array() + 3.0 * e11_gauss.array()) - 0.5 * P.array()).matrix().transpose());
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
