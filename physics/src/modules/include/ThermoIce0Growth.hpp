/*!
 * @file ThermoIce0Growth.hpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOICE0GROWTH_HPP
#define THERMOICE0GROWTH_HPP

#include "include/IVerticalIceGrowth.hpp"

namespace Nextsim {

class ThermoIce0Growth : public IVerticalIceGrowth {
public:
    ThermoIce0Growth()
        : IVerticalIceGrowth()
    {
        ModelComponent::requestProtectedArray(ModelComponent::ProtectedArray::H_ICE, &oldHi);
    }
    virtual ~ThermoIce0Growth() = default;

    void update(const TimestepTime& tsTime) override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField deltaHi;
    HField snowMelt;
    HField topMelt;
    HField botMelt;
    pConstHField oldHi;

    bool doFlooding = true; // TODO: read from configuration
};

} /* namespace Nextsim */

#endif /* THERMOICE0GROWTH_HPP */
