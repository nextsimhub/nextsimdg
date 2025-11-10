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
    void calculateElement(size_t i, const TimestepTime& tst);

    AdvectedField tInternal;
    AdvectedField tBottom;
    HField snowMelt;
    HField topMelt;
    HField botMelt;
    ModelArrayAccessor<Protected::SW_IN> sw_inAccessor;
    ModelArrayAccessor<Shared::SUBLIM, RO> sublAccessor;

    static const double cVol;
    static bool doFlooding;
    static const double seaIceTf;
    static double kappa_s;

    void calculateTemps(
        double& tSurf, double& tMidt, double& tBotn, double& mSurf, const double cice, const double dQia_dt, const double hice, const double hsnow, const double penSw, const double qia, const double tf, double dt);
};

} /* namespace Nextsim */

#endif /* THERMOWINTON_HPP */
