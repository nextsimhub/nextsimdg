/*!
 * @file ThermoIce0.hpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_THERMOICE0_HPP_
#define SRC_INCLUDE_THERMOICE0_HPP_

#include "include/Configured.hpp"
#include "IThermodynamics.hpp"

namespace Nextsim {

class PrognosticData;
class PhysicsData;
class ExternalData;
class NextsimPhysics;

//! The implementation class for the NeXtSIM therm0 ice thermodynamics.
class ThermoIce0 : public IThermodynamics, public Configured<ThermoIce0> {
public:
    ThermoIce0() = default;
    virtual ~ThermoIce0() = default;

    void configure() override;
    enum {
        KS_KEY,
        FLOODING_KEY,
    };

    /*!
     * @brief Calculate the NeXtSIM thermo0 ice thermodynamics.
     *
     * @param prog PrognosticData for this element (constant)
     * @param exter ExternalData for this element (constant)
     * @param phys PhysicsData for this element
     * @param nsphys Nextsim physics implementation data for this element.
     */
    void calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys,
        NextsimPhysics& nsphys) override;

private:
    static double k_s;
    static bool doFlooding;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_THERMOICE0_HPP_ */
