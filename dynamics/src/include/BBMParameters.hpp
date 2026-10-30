/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MEBPARAMETERS_HPP
#define __MEBPARAMETERS_HPP

#include "DynamicsParameters.hpp"

namespace Nextsim {

// static FloatType compactionParamDefault = -20.; //!< Compation parameter
static FloatType nu0Default = 1. / 3.; //!< \param Poisson's ratio
static FloatType youngDefault = 5.96e8; //!< \param Young's modulus
static FloatType P0Default = 10e3; //! < \param Constant to parametrize Pmax
static FloatType lambda0Default = 1e7; //!< \param lambda
static int alphaDefault = 5;
static FloatType expPMaxDefault = 1.5; //! \param Power of ice thickness in the pressure coefficient
static FloatType muDefault = 0.7; //!< \param tan_phi (FloatType) Internal friction coefficient (mu)
static FloatType comprCapDefault
    = 1e10; //! \param compr_strength (FloatType) Maximum compressive strength [N/m2]
static FloatType cLabDefault = 2e6; //! \param C_lab (FloatType) Test [Pa]
// static const int nStepsDefault = 120; //!< Number of sub-steps

class BBMParameters : public DynamicsParameters {

public:
    BBMParameters() = default;

    FloatType compactionParam = compactionParamDefault;
    FloatType nu0 = nu0Default;
    FloatType young = youngDefault;
    FloatType P0 = P0Default;
    FloatType lambda0 = lambda0Default;
    int alpha = alphaDefault;
    FloatType expPMax = expPMaxDefault;
    FloatType mu = muDefault;
    FloatType comprCap = comprCapDefault;
    FloatType cLab = cLabDefault;
    int nSteps = nStepsDefault;

    FloatType c0 = 10e3; //! \param

    FloatType minDamage = std::is_same_v<FloatType, float> ? 1e-6 : 1e-12;
};

} /* namespace Nextsim */

#endif /* __MEB_HPP */
