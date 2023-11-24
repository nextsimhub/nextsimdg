/*!
 * @file FluxConfiguredAtmosphere.hpp
 *
 * @date Sep 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FLUXCONFIGUREDATMOSPHERE_HPP
#define FLUXCONFIGUREDATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

#include "include/Configured.hpp"

namespace Nextsim {

//! A class to provide constant atmospheric forcings that can be configured at run time.
class FluxConfiguredAtmosphere : public IAtmosphereBoundary,
                                 public Configured<FluxConfiguredAtmosphere> {
public:
    FluxConfiguredAtmosphere() = default;
    ~FluxConfiguredAtmosphere() = default;

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
    std::string getName() const override { return "FluxConfiguredAtmosphere"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

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

#endif /* FLUXCONFIGUREDATMOSPHERE_HPP */
