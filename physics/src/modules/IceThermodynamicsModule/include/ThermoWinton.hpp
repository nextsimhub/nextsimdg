/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOWINTON_HPP
#define THERMOWINTON_HPP

#include "include/Configured.hpp"
#include "include/IIceThermodynamics.hpp"

namespace Nextsim {

//! A class implementing IIceThermodynamics as the Winton thermodynamics model.
class ThermoWinton : public IIceThermodynamics, public Configured<ThermoWinton> {
public:
    static const size_t nLevels;

    ThermoWinton();
    virtual ~ThermoWinton() = default;

    enum {
        KS_KEY,
        I0_KEY,
        FLOODING_KEY,
    };
    void configure() override;
    ConfigMap getConfiguration() const override;

    ModelState getStateDiagnostic() const override;
    ModelState getStatePrognostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tst) override;

    static const std::string tInteriorName;
    static const std::string tBottomName;

private:
    // local namespace to prevent conflicts with other thermodynamics implementations
    struct Private {
        static inline constexpr TextTag T_INTERNAL = "T_INTERNAL";
        static inline constexpr TextTag T_BOTTOM = "T_BOTTOM";
        static inline constexpr TextTag T_SNOW_MELT = "T_SNOW_MELT";
        static inline constexpr TextTag T_TOP_MELT = "T_TOP_MELT";
        static inline constexpr TextTag T_BOT_MELT = "T_BOT_MELT";
    };
    // private owned
    ModelArrayAccessor<Private::T_INTERNAL, RW> tInternalAccessor;
    ModelArrayAccessor<Private::T_BOTTOM, RW> tBottomAccessor;
    ModelArrayAccessor<Private::T_SNOW_MELT, RW> snowMeltAccessor;
    ModelArrayAccessor<Private::T_TOP_MELT, RW> topMeltAccessor;
    ModelArrayAccessor<Private::T_BOT_MELT, RW> botMeltAccessor;

    ModelArrayAccessor<Protected::SW_IN> sw_inAccessor;
    ModelArrayAccessor<Shared::SUBLIM, RO> sublAccessor;

    static const double cVol;
    static bool doFlooding;
    static const double seaIceTf;
    static double kappa_s;
};

} /* namespace Nextsim */

#endif /* THERMOWINTON_HPP */
