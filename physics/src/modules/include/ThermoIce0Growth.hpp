/*!
 * @file ThermoIce0Growth.hpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_MODULES_INCLUDE_THERMOICE0GROWTH_HPP_
#define PHYSICS_SRC_MODULES_INCLUDE_THERMOICE0GROWTH_HPP_

#include "include/IVerticalIceGrowth.hpp"

namespace Nextsim {

class ThermoIce0Growth : public IVerticalIceGrowth {
public:
    ThermoIce0Growth()
        : IVerticalIceGrowth()
    {
        ModelModule::requestProtectedArray(ModelModule::ProtectedArray::H_ICE, &oldHi);
    }
    virtual ~ThermoIce0Growth() = default;

    void update(const TimestepTime& tsTime) override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField deltaHi;
    HField snowMelt;
    HField topMelt;
    HField botMelt;
    const HField* oldHi;

    bool doFlooding = true; // TODO: read from configuration
};

} /* namespace Nextsim */

#endif /* PHYSICS_SRC_MODULES_INCLUDE_THERMOICE0GROWTH_HPP_ */
