/*!
 * @file ThermoIce0Growth.hpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOICE0GROWTH_HPP
#define THERMOICE0GROWTH_HPP

#include "include/IVerticalIceGrowth.hpp"

#include "include/Configured.hpp"
#include "include/ModelArrayRef.hpp"
namespace Nextsim {

class ThermoIce0Growth : public IVerticalIceGrowth, public Configured<ThermoIce0Growth> {
public:
    ThermoIce0Growth()
        : IVerticalIceGrowth()
    {
    }
    virtual ~ThermoIce0Growth() = default;

    enum {
        KS_KEY,
    };
    void configure() override;

    void update(const TimestepTime& tsTime) override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField snowMelt;
    HField topMelt;
    HField botMelt;
    HField qic;
    ModelArrayRef<ProtectedArray::H_ICE> oldHi;

    static double k_s;
    static const double freezingPointIce;

    bool doFlooding = true; // TODO: read from configuration
};

} /* namespace Nextsim */

#endif /* THERMOICE0GROWTH_HPP */
