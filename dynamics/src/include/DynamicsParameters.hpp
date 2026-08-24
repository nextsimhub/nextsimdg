/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 *
 */

#ifndef DYNAMICSPARAMETERS_HPP
#define DYNAMICSPARAMETERS_HPP

#include "NextsimDynamics.hpp"

namespace Nextsim {

static const FloatType rhoIceDefault = 900.; //!< Sea ice density
static const FloatType rhoAtmDefault = 1.3; //!< Air density
static const FloatType rhoOceanDefault = 1026.; //!< Ocean density
static const FloatType CAtmDefault = 1.2e-3; //!< Air drag coefficient
static const FloatType COceanDefault = 5.5e-3; //!< Ocean drag coefficient
static const FloatType fcDefault = 1.45842e-4; //!< Coriolis
static const FloatType oceanTurningAngleDefault = 25.; //!< Ocean turning angle

static const int nStepsDefault = 120; //!< number of sub-cycling steps
static const FloatType compactionParamDefault
    = -20.; //!< Compation parameter: Hibler's C in exp(-C(1-a))

class DynamicsParameters {

public:
    DynamicsParameters() = default;

    FloatType rhoIce = rhoIceDefault;
    FloatType rhoAtm = rhoAtmDefault;
    FloatType rhoOcean = rhoOceanDefault;

    FloatType CAtm = CAtmDefault;
    FloatType COcean = COceanDefault;

    FloatType fc = fcDefault;

    FloatType oceanTurningAngle = oceanTurningAngleDefault;
};
}

#endif /* DYNAMICSPARAMETERS_HPP */
