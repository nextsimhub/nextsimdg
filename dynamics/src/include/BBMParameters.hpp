/*!
 * @file BBMParameters.hpp
 * @date 19 Nov 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MEBPARAMETERS_HPP
#define __MEBPARAMETERS_HPP

#include "DynamicsParameters.hpp"

namespace Nextsim {

// static double compactionParamDefault = -20.; //!< Compation parameter
static double nu0Default = 1. / 3.; //!< \param Poisson's ratio
static double youngDefault = 5.96e8; //!< \param Young's modulus
static double P0Default = 10e3; //! < \param Constant to parametrize Pmax
static double lambda0Default = 1e7; //!< \param lambda
static int alphaDefault = 5;
static double expPMaxDefault = 1.5; //! \param Power of ice thickness in the pressure coefficient
static double muDefault = 0.7; //!< \param tan_phi (double) Internal friction coefficient (mu)
static double comprCapDefault
    = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]
static double cLabDefault = 2e6; //! \param C_lab (double) Test [Pa]
// static const int nStepsDefault = 120; //!< Number of sub-steps

class BBMParameters : public DynamicsParameters {

public:
    BBMParameters() = default;

    double compactionParam = compactionParamDefault;
    double nu0 = nu0Default;
    double young = youngDefault;
    double P0 = P0Default;
    double lambda0 = lambda0Default;
    int alpha = alphaDefault;
    double expPMax = expPMaxDefault;
    double mu = muDefault;
    double comprCap = comprCapDefault;
    double cLab = cLabDefault;
    int nSteps = nStepsDefault;

    double c0 = 10e3; //! \param
};

} /* namespace Nextsim */

#endif /* __MEB_HPP */
