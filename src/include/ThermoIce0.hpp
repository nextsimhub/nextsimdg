/*!
 * @file ThermoIce0.hpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_THERMOICE0_HPP_
#define SRC_INCLUDE_THERMOICE0_HPP_

#include "IThermodynamics.hpp"

namespace Nextsim {

class PrognosticData;
class PhysicsData;
class ExternalData;
class NextsimPhysics;

class ThermoIce0 : public IThermodynamics {
public:
    ThermoIce0() = default;
    virtual ~ThermoIce0() = default;

    void calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys,
        NextsimPhysics& nsphys);
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_THERMOICE0_HPP_ */
