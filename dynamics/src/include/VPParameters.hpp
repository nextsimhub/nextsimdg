/*!
 * @file VPParameters.hpp
 * @date 18 Jan 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __VPPARAMETERS_HPP
#define __VPPARAMETERS_HPP

#include "CommonDynamicsParameters.hpp"

namespace Nextsim {

class VPParameters : public CommonDynamicsParameters {
public:
    double Pstar; //!< Ice strength
    double DeltaMin; //!< Viscous regime

    VPParameters()
    {
        Pstar = 27500.0; //!< Ice strength
        DeltaMin = 2.e-9; //!< Viscous regime
    }
};

} /* namespace Nextsim */

#endif /* __VP_HPP */
