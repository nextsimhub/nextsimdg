/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOICE0HPP
#define THERMOICE0HPP

#include "include/Configured.hpp"
#include "include/IIceThermodynamics.hpp"

namespace Nextsim {

//! A class implementing IIceThermodynamics as the ThermoIce0 model.
class ThermoIce0 : public IIceThermodynamics, public Configured<ThermoIce0> {
public:
    ThermoIce0();
    virtual ~ThermoIce0() = default;

    enum {
        KS_KEY,
    };
    void configure() override;
    ConfigMap getConfiguration() const override;

    ModelState getStateDiagnostic() const override;
    ModelState getStatePrognostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tsTime) override;

private:
    // local namespace to prevent conflicts with other thermodynamics implementations
    struct Private {
        static inline constexpr TextTag T_SNOW_MELT = "T_SNOW_MELT";
        static inline constexpr TextTag T_TOP_MELT = "T_TOP_MELT";
        static inline constexpr TextTag T_BOT_MELT = "T_BOT_MELT";
    };
    // private owned
    ModelArrayAccessor<Private::T_SNOW_MELT, RW> snowMeltAccessor;
    ModelArrayAccessor<Private::T_TOP_MELT, RW> topMeltAccessor;
    ModelArrayAccessor<Private::T_BOT_MELT, RW> botMeltAccessor;

    ModelArrayAccessor<Shared::Q_IC, RW> qicAccessor;

    static const double freezingPointIce;
    static double kappa_s;

    bool doFlooding = true; // TODO: read from configuration
};

} /* namespace Nextsim */

#endif /* THERMOICE0HPP */
