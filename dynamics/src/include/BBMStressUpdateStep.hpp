/*!
 * @file BBMStressUpdateStep.hpp
 *
 * @date 13 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BBMSTRESSUPDATESTEP_HPP
#define BBMSTRESSUPDATESTEP_HPP

#include "include/StressUpdateStep.hpp"

#include "include/BBMParameters.hpp"
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
    ~BBMStressUpdateStep() override = default;
    void stressUpdateHighOrder(const DynamicsParameters& dParams, const ParametricMesh& smesh,
        SymmetricTensorVector& stress, const SymmetricTensorVector& strain,
        const DGVector<DGadvection>& h, const DGVector<DGadvection>& a,
        const double deltaT) override
    {
        assert(pmap);

        // Unwrap references
        DGVector<DGstress>& s11 = stress[I11];
        DGVector<DGstress>& s12 = stress[I12];
        DGVector<DGstress>& s22 = stress[I22];

        DGVector<DGstress>& e11 = strain[I11];
        DGVector<DGstress>& e12 = strain[I12];
        DGVector<DGstress>& e22 = strain[I22];

        const auto& params = reinterpret_cast<const BBMParameters&>(dParams);
        // Number of Gauss points
        const size_t nGauss1D = GAUSSPOINTS1D(DGstress);
        const size_t nGauss2D = GAUSSPOINTS(DGstress);

        const double sqrtNuRhoI = std::sqrt(2. * (1. + params.nu0) * params.rhoIce);

//! Stress and Damage Update
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nelements; ++i) {

            if (!smesh.landmask[i])
                continue;

            //! Evaluate values in Gauss points (3 point Gauss rule in 2d => 9 points)
            const LocalEdgeVector<nGauss2D> hGauss
                = (h.row(i) * PSI<DGadvection, nGauss1D>).array().max(0.0);
            const LocalEdgeVector<nGauss2D> aGauss
                = (a.row(i) * PSI<DGadvection, nGauss1D>).array().max(0.0).min(1.0);

            LocalEdgeVector<nGauss2D> dGauss
                = (p_d->row(i) * PSI<DGadvection, nGauss1D>).array().max(1e-12).min(1.0);

            const LocalEdgeVector<nGauss2D> e11Gauss = e11.row(i) * PSI<DGstress, nGauss1D>;
            const LocalEdgeVector<nGauss2D> e12Gauss = e12.row(i) * PSI<DGstress, nGauss1D>;
            const LocalEdgeVector<nGauss2D> e22Gauss = e22.row(i) * PSI<DGstress, nGauss1D>;

            LocalEdgeVector<nGauss2D> s11Gauss = s11.row(i) * PSI<DGstress, nGauss1D>;
            LocalEdgeVector<nGauss2D> s12Gauss = s12.row(i) * PSI<DGstress, nGauss1D>;
            LocalEdgeVector<nGauss2D> s22Gauss = s22.row(i) * PSI<DGstress, nGauss1D>;

            //! Current normal stress for the evaluation of tildeP (Eqn. 1)
            LocalEdgeVector<nGauss2D> sigmaN = 0.5 * (s11Gauss + s22Gauss).array();

            //! exp(-C(1-A))
            const LocalEdgeVector<nGauss2D> expC
                = (params.compactionParam * (1. - aGauss.array())).exp();

            // Eqn. 25
            const LocalEdgeVector<nGauss2D> viscousTime
                = params.lambda0 * (dGauss.array() * expC.array()).pow(params.alpha - 1);

            //! BBM  Computing tildeP according to (Eqn. 7b and Eqn. 8)
            // (Eqn. 8)
            // NB! We need to add 1 to the exponent because stress is thickness integrated
            const LocalEdgeVector<nGauss2D> Pmax
                = params.P0 * hGauss.array().pow(params.expPMax + 1.) * expC.array();

            // (Eqn. 7b) Prepare tildeP
            // tildeP must be capped at 1 to get an elastic response
            // (Eqn. 7b) Select case based on sigma_n
            const LocalEdgeVector<nGauss2D> tildeP
                = (sigmaN.array() < 0.).select((-Pmax.array() / sigmaN.array()).min(1.), 0.);

            // multiplicator
            const LocalEdgeVector<nGauss2D> multiplicator
                = viscousTime.array() / (viscousTime.array() + (1. - tildeP.array()) * deltaT);

            //! Eqn. 9
            const LocalEdgeVector<nGauss2D> elasticity
                = params.young * dGauss.array() * expC.array();

            // Eqn. 12: first factor on RHS
            /* Stiffness matrix
             * / (K:e)11 \       1     /  1  nu    0  \ / e11 \
             * | (K:e)22 |  =  ------- | nu   1    0  | | e22 |
             * \ (K:e)12 /    1 - nu^2 \  0   0  1-nu / \ e12 /
             */

            // NB! We multiply with h because stress is thickness integrated
            const LocalEdgeVector<nGauss2D> KFactor
                = deltaT * hGauss.array() * elasticity.array() / (1. - params.nu0 * params.nu0);

            s11Gauss.array() += KFactor.array() * (e11Gauss + params.nu0 * e22Gauss).array();
            s22Gauss.array() += KFactor.array() * (params.nu0 * e11Gauss + e22Gauss).array();
            s12Gauss.array() += KFactor.array() * e12Gauss.array() * (1. - params.nu0);

            //! Implicit part of RHS (Eqn. 33)
            s11Gauss.array() *= multiplicator.array();
            s22Gauss.array() *= multiplicator.array();
            s12Gauss.array() *= multiplicator.array();

            sigmaN = 0.5 * (s11Gauss + s22Gauss).array();
            const LocalEdgeVector<nGauss2D> tau
                = (0.25 * (s11Gauss - s22Gauss).array().square() + s12Gauss.array().square())
                      .sqrt();

            //! Eqn. 30
            const LocalEdgeVector<nGauss2D> comprStrength
                = params.comprCap * cohesion[i] * hGauss.array();

            // Compressive and Mohr-Coulomb failure using Messrs Plante & Tremblay's formulation
            LocalEdgeVector<nGauss2D> dCrit
                = (sigmaN.array() < -comprStrength.array())
                      .select(-comprStrength.array() / sigmaN.array(),
                          cohesion[i] * hGauss.array()
                              / (tau.array() + params.mu * sigmaN.array()));

            // Only damage when we're outside, i.e. dCrit < 1
            // tau - mu*sigmaN < 0 is always inside but gives dCrit < 0
            dCrit = (dCrit.array() < 0.).select(1.0, dCrit.array().min(1.));

            // Eqn. 29 (but we calculate the reciprocal directly)
            const LocalEdgeVector<nGauss2D> rtd
                = elasticity.array().sqrt() / (smesh.h(i) * sqrtNuRhoI);

            // Update damage
            dGauss.array() -= dGauss.array() * (1. - dCrit.array()) * deltaT * rtd.array();

            // Relax stress in Gassus points
            s11Gauss.array() -= s11Gauss.array() * (1. - dCrit.array()) * deltaT * rtd.array();
            s12Gauss.array() -= s12Gauss.array() * (1. - dCrit.array()) * deltaT * rtd.array();
            s22Gauss.array() -= s22Gauss.array() * (1. - dCrit.array()) * deltaT * rtd.array();

            // INTEGRATION OF STRESS AND DAMAGE
            s11.row(i) = pmap->iMJwPSI[i] * s11Gauss.matrix().transpose();
            s12.row(i) = pmap->iMJwPSI[i] * s12Gauss.matrix().transpose();
            s22.row(i) = pmap->iMJwPSI[i] * s22Gauss.matrix().transpose();
            p_d->row(i) = pmap->iMJwPSI_dam[i] * dGauss.matrix().transpose();
        }
    }

    void setDamage(DGVector<DGadvection>& dIn) { p_d = &dIn; }
    void setPMap(ParametricMomentumMap<CG, DGadvection>* pmapIn) { pmap = pmapIn; }
    void setCohesion(const BBMParameters& params, const ParametricMesh& smesh)
    {
        cohesion.resize(smesh.nelements);
        for (size_t i = 0; i < smesh.nelements; ++i) {
            //! Eqn. 22
            const double scaleCoef
                = params.scaleCohesion ? std::sqrt(params.referenceScaleC / smesh.h(i)) : 1.;
            cohesion[i] = params.cohesion * scaleCoef;
        }
    }

protected:
    ParametricMomentumMap<CG, DGadvection>* pmap;
    DGVector<DGadvection>* p_d;
    std::vector<double> cohesion;
};
}

#endif /* BBMSTRESSUPDATESTEP_ */
