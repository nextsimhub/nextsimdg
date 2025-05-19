/*!
 * @file BBMParameters.hpp
 *
 * @date 19 May 2025
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
static double comprCapDefault = 1e4; //! \param comprCap (double) Maximum compressive strength
                                     //! multiplier (cap as fraction of cohesion).
static double cohesionDefault
    = 2e6; //! \param C_lab (double) Cohesion (at the lab scale if scaling is turned on) [Pa]
static bool scaleCohesionDefault = true; //! Scale the cohesion as sqrt of the grid size
static double referenceScaleCDefault = 0.1; //! \param Reference scale for cohesion

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
    double cohesion = cohesionDefault;
    bool scaleCohesion = scaleCohesionDefault;
    double referenceScaleC = referenceScaleCDefault;

    double c0 = 10e3; //! \param
};

} /* namespace Nextsim */

#endif /* __MEB_HPP */
