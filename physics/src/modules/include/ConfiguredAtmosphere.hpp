/*!
 * @file ConfiguredAtmosphere.hpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDATMOSPHERE_HPP
#define CONFIGUREDATMOSPHERE_HPP

#include "AtmosphereState.hpp"

namespace Nextsim {

class ConfiguredAtmosphere : public AtmosphereState, public Configured<ConfiguredAtmosphere> {
public:
    ConfiguredAtmosphere() = default;
    ~ConfiguredAtmosphere() = default;

    enum {
        TAIR_KEY,
        TDEW_KEY,
        PAIR_KEY,
        RMIX_KEY,
        SWIN_KEY,
        LWIN_KEY,
        SNOW_KEY,
        WIND_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ConfiguredAtmosphere"; }

    void configure() override;

protected:
    //! Performs the implementation specific updates. Does nothing.
    void updateSpecial(const TimestepTime&) override { }

private:
    static double tair0;
    static double tdew0;
    static double pair0;
    static double rmix0;
    static double sw_in0;
    static double lw_in0;
    static double snowfall0;
    static double windspeed0;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDATMOSPHERE_HPP */
