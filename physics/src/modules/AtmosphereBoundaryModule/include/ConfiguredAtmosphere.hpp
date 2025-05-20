/*!
 * @file ConfiguredAtmosphere.hpp
 *
 * @date 30 Apr 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDATMOSPHERE_HPP
#define CONFIGUREDATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

#include "include/IFluxCalculation.hpp"

namespace Nextsim {

//! A class to provide constant atmospheric forcings that can be configured at run time.
class ConfiguredAtmosphere : public IAtmosphereBoundary, public Configured<ConfiguredAtmosphere> {
public:
    ConfiguredAtmosphere();
    ~ConfiguredAtmosphere() = default;

    enum {
        TAIR_KEY,
        TDEW_KEY,
        PAIR_KEY,
        SW_KEY,
        LW_KEY,
        SNOW_KEY,
        RAIN_KEY,
        UWIND_KEY,
        VWIND_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ConfiguredAtmosphere"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;
    ConfigMap getConfiguration() const override;

    //! Calculates the fluxes from the given values
    void update(const TimestepTime&) override;

private:
    static double tair0;
    static double tdew0;
    static double pair0;
    static double sw0;
    static double lw0;
    static double snowfall0;
    static double rain0;
    static double uWind0;
    static double vWind0;

    HField tair;
    HField tdew;
    HField pair;
    HField sw_in;
    HField lw_in;
    HField wind;

    IFluxCalculation* fluxImpl;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDATMOSPHERE_HPP */
