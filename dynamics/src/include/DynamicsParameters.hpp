/*!
 * @file DynamicsParameters.hpp
 *
 * @date 19 Nov 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 */

#ifndef DYNAMICSPARAMETERS_HPP
#define DYNAMICSPARAMETERS_HPP

namespace Nextsim {

static const double rhoIceDefault = 900.; //!< Sea ice density
static const double rhoAtmDefault = 1.3; //!< Air density
static const double rhoOceanDefault = 1026.; //!< Ocean density
static const double CAtmDefault = 1.2e-3; //!< Air drag coefficient
static const double COceanDefault = 5.5e-3; //!< Ocean drag coefficient
static const double fcDefault = 1.45842e-4; //!< Coriolis
static const double oceanTurningAngleDefault = 25.; //!< Ocean turning angle

static const int nStepsDefault = 120; //!< number of sub-cycling steps
static const double compactionParamDefault
    = -20.; //!< Compation parameter: Hibler's C in exp(-C(1-a))

class DynamicsParameters {

public:
    DynamicsParameters() = default;

    double rhoIce = rhoIceDefault;
    double rhoAtm = rhoAtmDefault;
    double rhoOcean = rhoOceanDefault;

    double CAtm = CAtmDefault;
    double COcean = COceanDefault;

    double fc = fcDefault;

    double oceanTurningAngle = oceanTurningAngleDefault;
};
}

#endif /* DYNAMICSPARAMETERS_HPP */
