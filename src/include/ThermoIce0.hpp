/*!
 * @file ThermoIce0.hpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_THERMOICE0_HPP_
#define SRC_INCLUDE_THERMOICE0_HPP_

#include "Configured.hpp"
#include "IThermodynamics.hpp"

namespace Nextsim {

class PrognosticData;
class PhysicsData;
class ExternalData;
class NextsimPhysics;

class ThermoIce0 : public IThermodynamics, public Configured<ThermoIce0> {
public:
    ThermoIce0() = default;
    virtual ~ThermoIce0() = default;

    void configure() override;
    enum {
        KS_KEY,
        FLOODING_KEY,
    };

    void calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys,
        NextsimPhysics& nsphys);

private:
    static double k_s;
    static bool doFlooding;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_THERMOICE0_HPP_ */
