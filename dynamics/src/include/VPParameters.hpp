/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __VPPARAMETERS_HPP
#define __VPPARAMETERS_HPP

#include "DynamicsParameters.hpp"

namespace Nextsim {

// static const int nStepsDefault = 100; //!< number of sub-cycling steps
// static const FloatType compactionParamDefault = -20.; //!< Compation parameter: Hibler's C in
// exp(-C(1-a))
static const FloatType pStarDefault = 27.5e3; //!< Ice strength
static const FloatType deltaMinDefault = 2e-9; //!< Viscous regime

class VPParameters : public DynamicsParameters {

public:
    VPParameters() = default;

    FloatType compactionParam = compactionParamDefault;
    FloatType pStar = pStarDefault;
    FloatType deltaMin = deltaMinDefault;
    int nSteps = nStepsDefault;
};

} /* namespace Nextsim */

#endif /* __VP_HPP */
