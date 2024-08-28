/*!
 * @file BBMStressUpdateStep.hpp
 *
 * @date Mar 1, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMSTRESSUPDATESTEP_HPP
#define BBMSTRESSUPDATESTEP_HPP

#include "include/StressUpdateStep.hpp"

#include "include/MEBParameters.hpp"
#include "include/ParametricMap.hpp"
#include "include/codeGenerationDGinGauss.hpp"

namespace Nextsim {

template <int DGadvection, int DGstress, int CG>
class BBMStressUpdateStep : public StressUpdateStep<DGadvection, DGstress> {

public:
    typedef std::array<std::reference_wrapper<DGVector<DGstress>>, N_TENSOR_ELEMENTS>
        SymmetricTensorVector;
    BBMStressUpdateStep()
        : pmap(nullptr)
        , p_d(nullptr)
    {
    }
    ~BBMStressUpdateStep() = default;
    void stressUpdateHighOrder(const DynamicsParameters& dParams, const ParametricMesh& smesh,
        SymmetricTensorVector& stress, const SymmetricTensorVector& strain,
        const DGVector<DGadvection>& h, const DGVector<DGadvection>& a,
        const double deltaT) override
    {
        // Number of Gauss points
        static constexpr size_t NGP
            = (((DGstress == 8) || (DGstress == 6)) ? 3 : (DGstress == 3 ? 2 : -1));
        using EdgeVec = Eigen::Matrix<FloatType, 1, NGP * NGP>;

        // Unwrap references
        DGVector<DGstress>& s11 = stress[I11];
        DGVector<DGstress>& s12 = stress[I12];
        DGVector<DGstress>& s22 = stress[I22];

        DGVector<DGstress>& e11 = strain[I11];
        DGVector<DGstress>& e12 = strain[I12];
        DGVector<DGstress>& e22 = strain[I22];

        const MEBParameters& params = reinterpret_cast<const MEBParameters&>(dParams);

//! Stress and Damage Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            //! Evaluate values in Gauss points (3 point Gauss rule in 2d => 9 points)
            const EdgeVec hGauss = (h.row(i) * PSI<DGadvection, NGP>).array().max(0.0).matrix();
            const EdgeVec aGauss
                = (a.row(i) * PSI<DGadvection, NGP>).array().max(0.0).min(1.0).matrix();
            EdgeVec dGauss
                = (p_d->row(i) * PSI<DGadvection, NGP>).array().max(1e-12).min(1.0).matrix();

            const EdgeVec e11Gauss = e11.row(i) * PSI<DGstress, NGP>;
            const EdgeVec e12Gauss = e12.row(i) * PSI<DGstress, NGP>;
            const EdgeVec e22Gauss = e22.row(i) * PSI<DGstress, NGP>;

            EdgeVec s11Gauss = s11.row(i) * PSI<DGstress, NGP>;
            EdgeVec s12Gauss = s12.row(i) * PSI<DGstress, NGP>;
            EdgeVec s22Gauss = s22.row(i) * PSI<DGstress, NGP>;

            //! Current normal stress for the evaluation of tildeP (Eqn. 1)
            EdgeVec sigma_n = 0.5 * (s11Gauss.array() + s22Gauss.array());

            //! exp(-C(1-A))
            const EdgeVec expC = (params.compaction_param * (1.0 - aGauss.array())).exp().array();

            // Eqn. 25
            const EdgeVec powalphaexpC
                = (dGauss.array() * expC.array()).pow(params.exponent_relaxation_sigma - 1);
            const EdgeVec time_viscous = params.undamaged_time_relaxation_sigma * powalphaexpC;

            //! BBM  Computing tildeP according to (Eqn. 7b and Eqn. 8)
            // (Eqn. 8)
            const EdgeVec Pmax
                = params.P0 * hGauss.array().pow(params.exponent_compression_factor) * expC.array();

            // (Eqn. 7b) Prepare tildeP
            // tildeP must be capped at 1 to get an elastic response
            // (Eqn. 7b) Select case based on sigma_n
            const EdgeVec tildeP
                = (sigma_n.array() < 0.0)
                      .select((-Pmax.array() / sigma_n.array()).min(1.0).matrix(), 0.);

            // multiplicator
            const EdgeVec multiplicator
                = time_viscous.array() / (time_viscous.array() + (1. - tildeP.array()) * deltaT);

            //! Eqn. 9
            const EdgeVec elasticity
                = hGauss.array() * params.young * dGauss.array() * expC.array();

            // Eqn. 12: first factor on RHS
            /* Stiffness matrix
             * / (K:e)11 \       1     /  1  nu    0  \ / e11 \
             * | (K:e)22 |  =  ------- | nu   1    0  | | e22 |
             * \ (K:e)12 /    1 - nu^2 \  0   0  1-nu / \ e12 /
             */

            const EdgeVec Dunit_factor
                = deltaT * elasticity.array() / (1. - (params.nu0 * params.nu0));

            s11Gauss.array()
                += Dunit_factor.array() * (e11Gauss.array() + params.nu0 * e22Gauss.array());
            s22Gauss.array()
                += Dunit_factor.array() * (params.nu0 * e11Gauss.array() + e22Gauss.array());
            s12Gauss.array() += Dunit_factor.array() * e12Gauss.array() * (1. - params.nu0);

            //! Implicit part of RHS (Eqn. 33)
            s11Gauss.array() *= multiplicator.array();
            s22Gauss.array() *= multiplicator.array();
            s12Gauss.array() *= multiplicator.array();

            sigma_n = 0.5 * (s11Gauss.array() + s22Gauss.array());
            const EdgeVec tau = (0.25 * (s11Gauss.array() - s22Gauss.array()).square()
                + s12Gauss.array().square())
                                    .sqrt();

            const FloatType scale_coef = std::sqrt(0.1 / smesh.h(i));

            //! Eqn. 22
            const EdgeVec cohesion = params.C_lab * scale_coef * hGauss.array();
            //! Eqn. 30
            const EdgeVec compr_strength = params.compr_strength * scale_coef * hGauss.array();

            // Mohr-Coulomb failure using Mssrs. Plante & Tremblay's formulation
            // sigma_s + tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
            EdgeVec dcrit
                = (tau.array() + params.tan_phi * sigma_n.array() > 0.)
                      .select(
                          cohesion.array() / (tau.array() + params.tan_phi * sigma_n.array()), 1.);

            // Compressive failure using Mssrs. Plante & Tremblay's formulation
            dcrit = (sigma_n.array() < -compr_strength.array())
                        .select(-compr_strength.array() / sigma_n.array(), dcrit);

            // Only damage when we're outside
            dcrit = dcrit.array().min(1.0);

            // Eqn. 29
            const EdgeVec td = smesh.h(i) * std::sqrt(2. * (1. + params.nu0) * params.rho_ice)
                / elasticity.array().sqrt();

            // Update damage
            dGauss.array() -= dGauss.array() * (1. - dcrit.array()) * deltaT / td.array();

            // Relax stress in Gassus points
            s11Gauss.array() -= s11Gauss.array() * (1. - dcrit.array()) * deltaT / td.array();
            s12Gauss.array() -= s12Gauss.array() * (1. - dcrit.array()) * deltaT / td.array();
            s22Gauss.array() -= s22Gauss.array() * (1. - dcrit.array()) * deltaT / td.array();

            // INTEGRATION OF STRESS AND DAMAGE
            // get the inverse of the mass matrix scaled with the test-functions in the gauss
            // points, with the gauss weights and with J. This is a 8 x 9 matrix
            const auto& iMJwPSI = pmap->iMJwPSI[i];
            s11.row(i) = iMJwPSI * s11Gauss.matrix().transpose();
            s12.row(i) = iMJwPSI * s12Gauss.matrix().transpose();
            s22.row(i) = iMJwPSI * s22Gauss.matrix().transpose();

            p_d->row(i) = pmap->iMJwPSIAdvect[i] * dGauss.matrix().transpose();
        }
    }

    void setDamage(DGVector<DGadvection>& dIn) { p_d = &dIn; }
    void setPMap(ParametricMomentumMap<CG>* pmapIn) { pmap = pmapIn; }

protected:
    ParametricMomentumMap<CG>* pmap;
    DGVector<DGadvection>* p_d;
};
}

#endif /* BBMSTRESSUPDATESTEP_ */
