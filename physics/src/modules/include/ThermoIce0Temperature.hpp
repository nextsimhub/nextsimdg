/*!
 * @file ThermoIce0Temperature.hpp
 *
 * @date Apr 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOICE0TEMPERATURE_HPP
#define THERMOICE0TEMPERATURE_HPP

#include "include/Configured.hpp"
#include "include/IIceTemperature.hpp"
#include "include/constants.hpp"

namespace Nextsim {

class ThermoIce0Temperature : public IIceTemperature, public Configured<ThermoIce0Temperature> {
public:
    void update(const TimestepTime& tst) override
    {
        overElements(std::bind(&ThermoIce0Temperature::calculateElement, this,
                         std::placeholders::_1, std::placeholders::_2),
            tst);
    }

    void configure() override;
    enum {
        KS_KEY,
    };

private:
    static double k_s;
    static const double freezingPointIce;

    void calculateElement(size_t i, const TimestepTime& tst);
};

} /* namespace Nextsim */

#endif /* THERMOICE0TEMPERATURE_HPP */
