/*!
 * @file ConfiguredAtmosphere.hpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDATMOSPHERE_HPP
#define CONFIGUREDATMOSPHERE_HPP

#include "IAtmosphereBoundary.hpp"

namespace Nextsim {

//! A class to provide constant atmospheric forcings that can be configured at run time.
class ConfiguredAtmosphere : public IAtmosphereBoundary, public Configured<ConfiguredAtmosphere> {
public:
    ConfiguredAtmosphere() = default;
    ~ConfiguredAtmosphere() = default;

    enum {
        QIA_KEY,
        DQIA_DT_KEY,
        QOW_KEY,
        SUBL_KEY,
        SNOW_KEY,
        RAIN_KEY,
        EVAP_KEY,
        WINDU_KEY,
        WINDV_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ConfiguredAtmosphere"; }

    void configure() override;

protected:
    //! Performs the implementation specific updates. Does nothing.
    void update(const TimestepTime&) override { }

private:
    static double qia0;
    static double dqia_dt0;
    static double qow0;
    static double subl0;
    static double snowfall0;
    static double rain0;
    static double evap0;
    static double u0;
    static double v0;

};

} /* namespace Nextsim */

#endif /* CONFIGUREDATMOSPHERE_HPP */
