/*!
 * @file ThermoIce0.hpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOICE0HPP
#define THERMOICE0HPP

#include "include/Configured.hpp"
#include "include/IIceThermodynamics.hpp"
#include "include/ModelArrayRef.hpp"
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

    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tsTime) override;

    size_t getNZLevels() const override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField snowMelt;
    HField topMelt;
    HField botMelt;
    HField qic;
    ModelArrayRef<Protected::HTRUE_ICE> oldHi;

    static const double freezingPointIce;
    static double kappa_s;

    bool doFlooding = true; // TODO: read from configuration

    static const size_t nZLevels;
};

} /* namespace Nextsim */

#endif /* THERMOICE0HPP */
