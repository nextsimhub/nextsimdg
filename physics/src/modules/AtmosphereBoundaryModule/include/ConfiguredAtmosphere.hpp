/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
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
    static FloatType tair0;
    static FloatType tdew0;
    static FloatType pair0;
    static FloatType sw0;
    static FloatType lw0;
    static FloatType snowfall0;
    static FloatType rain0;
    static FloatType uWind0;
    static FloatType vWind0;

    ModelArrayAccessor<Protected::T_AIR, RW> tairAccessor;
    ModelArrayAccessor<Protected::DEW_2M, RW> tdewAccessor;
    ModelArrayAccessor<Protected::P_AIR, RW> pairAccessor;
    ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;
    ModelArrayAccessor<Protected::LW_IN, RW> lw_inAccessor;
    ModelArrayAccessor<Protected::WIND_SPEED, RW> windAccessor;

    IFluxCalculation* fluxImpl;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDATMOSPHERE_HPP */
