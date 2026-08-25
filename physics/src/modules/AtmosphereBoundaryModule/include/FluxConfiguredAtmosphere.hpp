/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
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
    ConfigMap getConfiguration() const override;

protected:
    //! Performs the implementation specific updates.
    void update(const TimestepTime&) override;

private:
    static FloatType qia0;
    static FloatType dqia_dt0;
    static FloatType qow0;
    static FloatType subl0;
    static FloatType snowfall0;
    static FloatType rain0;
    static FloatType evap0;
    static FloatType u0;
    static FloatType v0;
};

} /* namespace Nextsim */

#endif /* FLUXCONFIGUREDATMOSPHERE_HPP */
