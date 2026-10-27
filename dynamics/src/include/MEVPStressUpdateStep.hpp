/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MEVPSTRESSUPDATESTEP_HPP
#define MEVPSTRESSUPDATESTEP_HPP

#include "include/StressUpdateStep.hpp"

#include "include/ParametricMap.hpp"
#include "include/VPParameters.hpp"
#include "include/codeGenerationDGinGauss.hpp"

namespace Nextsim {

/*!
 * The mEVP stress update in all cells:
 *   S = (1-1/alpha) S + 1/alpha S_VP(E)
 * iMJwPSI[i] is M^-1 J w PSI of the stress space in the GAUSSPOINTS1D(DGstress)^2 Gauss
 * points of cell i. Used by the CG momentum (MEVPStressUpdateStep) and the LDG momentum.
 */
template <int DGadvection, int DGstress, class GaussMapMatrixVector>
void mEVPStressUpdate(const VPParameters& vpParams, const ParametricMesh& smesh,
    const GaussMapMatrixVector& iMJwPSI, DGVector<DGstress>& s11, DGVector<DGstress>& s12,
    DGVector<DGstress>& s22, const DGVector<DGstress>& e11, const DGVector<DGstress>& e12,
    const DGVector<DGstress>& e22, const DGVector<DGadvection>& h, const DGVector<DGadvection>& a,
    const FloatType alpha)
{
    const FloatType sqrDeltaMin = SQR(vpParams.deltaMin);

    // Number of Gauss points
    constexpr int nGauss = GAUSSPOINTS1D(DGstress);

    //! Stress Update
#pragma omp parallel for
    for (size_t i = 0; i < smesh.nelements; ++i) {
        if (!smesh.landmask[i])
            continue;

        // Here, one should check if it is enough to use a 2-point Gauss rule.
        // We're dealing with dG2, 3-point Gauss should be required.

        const LocalEdgeVector<nGauss * nGauss> h_gauss
            = (h.row(i) * PSI<DGadvection, nGauss>).array().max(0.0_ft).matrix();
        const LocalEdgeVector<nGauss * nGauss> a_gauss
            = (a.row(i) * PSI<DGadvection, nGauss>).array().max(0.0_ft).min(1.0_ft).matrix();

        const LocalEdgeVector<nGauss * nGauss> e11_gauss = e11.row(i) * PSI<DGstress, nGauss>;
        const LocalEdgeVector<nGauss * nGauss> e12_gauss = e12.row(i) * PSI<DGstress, nGauss>;
        const LocalEdgeVector<nGauss * nGauss> e22_gauss = e22.row(i) * PSI<DGstress, nGauss>;

        const LocalEdgeVector<nGauss * nGauss> DELTA
            = (sqrDeltaMin + 1.25_ft * (e11_gauss.array().square() + e22_gauss.array().square())
                + 1.50_ft * e11_gauss.array() * e22_gauss.array() + e12_gauss.array().square())
                  .sqrt()
                  .matrix();
        // FloatType DELTA = sqrt(SQR(vpparameters.deltaMin) + 1.25 * (SQR(E11(i, 0)) +
        // SQR(E22(i, 0)))
        //       + 1.50 * E11(i, 0) * E22(i, 0) + SQR(E12(i, 0)));
        //   assert(DELTA > 0);

        //   //! Ice strength
        //   FloatType P = vpparameters.pStar * H(i, 0) * exp(-20.0 * (1.0 - A(i, 0)));
        const LocalEdgeVector<nGauss * nGauss> P = (vpParams.pStar * h_gauss.array()
            * (vpParams.compactionParam * (1.0_ft - a_gauss.array())).exp())
                                                       .matrix();

        // replacement pressure P_r = P * DELTA / (DELTA + deltaMin), the viscosities keep P
        const LocalEdgeVector<nGauss * nGauss> Pr = vpParams.replacementPressure
            ? (P.array() * DELTA.array() / (DELTA.array() + vpParams.deltaMin)).matrix()
            : P;

        //   // S = S_old + 1/alpha (S(u)-S_old) = (1-1/alpha) S_old + 1/alpha S(u)
        s11.row(i) *= FloatType(1.0 - 1.0 / alpha);
        s12.row(i) *= FloatType(1.0 - 1.0 / alpha);
        s22.row(i) *= FloatType(1.0 - 1.0 / alpha);

        // const Eigen::Matrix<FloatType, 1, 9> J = ParametricTools::J<3>(smesh, i);
        // // get the inverse of the mass matrix scaled with the test-functions in the gauss
        // points,
        // // with the gauss weights and with J. This is a 8 x 9 matrix
        // const Eigen::Matrix<FloatType, 8, 9> imass_psi =
        // ParametricTools::massMatrix<8>(smesh, i).inverse()
        //     * (PSI<8,3>.array().rowwise() * (GAUSSWEIGHTS<3>.array() * J.array())).matrix();

        s11.row(i) += iMJwPSI[i]
            * (1.0_ft / alpha
                * (P.array() / 8.0_ft / DELTA.array()
                        * (5.0_ft * e11_gauss.array() + FloatType(3) * e22_gauss.array())
                    - 0.5_ft * Pr.array())
                    .matrix()
                    .transpose());

        //   S12.row(i) += 1.0 / alpha * (2. * eta * E12.row(i));
        // 2 eta = 2/4 * P / (2 Delta) = P / (4 Delta)
        s12.row(i) += iMJwPSI[i]
            * (1.0_ft / alpha
                * (P.array() / FloatType(4) / DELTA.array() * e12_gauss.array())
                    .matrix()
                    .transpose());

        //   S22.row(i)
        //       += 1.0 / alpha * (2. * eta * E22.row(i) + (zeta - eta) * (E11.row(i) +
        //       E22.row(i)));
        //   S22(i, 0) -= 1.0 / alpha * 0.5_ft * P;
        s22.row(i) += iMJwPSI[i]
            * (1.0_ft / alpha
                * (P.array() / 8.0_ft / DELTA.array()
                        * (5.0_ft * e22_gauss.array() + FloatType(3) * e11_gauss.array())
                    - 0.5_ft * Pr.array())
                    .matrix()
                    .transpose());
    }
}

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
        const DGVector<DGadvection>& h, const DGVector<DGadvection>& a,
        const FloatType deltaT) override
    {
        mEVPStressUpdate<DGadvection, DGstress>(reinterpret_cast<const VPParameters&>(params),
            smesh, pmap->iMJwPSI, stress[I11], stress[I12], stress[I22], strain[I11], strain[I12],
            strain[I22], h, a, alpha);
    }

    void setPMap(ParametricMomentumMap<CG, DGadvection>* pmapIn) { pmap = pmapIn; }

protected:
    ParametricMomentumMap<CG, DGadvection>* pmap;

private:
    //! MEVP parameters
    FloatType alpha = 1500.0;
    FloatType beta = 1500.0;
};

} /* namespace Nextsim */

#endif /* MEVPSTRESSUPDATESTEP_HPP */
