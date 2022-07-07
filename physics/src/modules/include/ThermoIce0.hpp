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

//! A class implementing IICeThermodynamics as the ThermoIce0 model.
class ThermoIce0 : public IIceThermodynamics, public Configured<ThermoIce0> {
public:
    ThermoIce0();
    virtual ~ThermoIce0() = default;

    enum {
        KS_KEY,
    };
    void configure() override;

    void setData(const ModelState&) override;
    void update(const TimestepTime& tsTime) override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField snowMelt;
    HField topMelt;
    HField botMelt;
    HField qic;
    ModelArrayRef<ProtectedArray::HTRUE_ICE> oldHi;

    static double k_s;
    static const double freezingPointIce;

    bool doFlooding = true; // TODO: read from configuration
};

} /* namespace Nextsim */

#endif /* THERMOICE0HPP */
