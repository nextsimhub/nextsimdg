/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/MU71Atmosphere.hpp"

#include "include/Finalizer.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/NextsimModule.hpp"

namespace Nextsim {

void MU71Atmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);

    iIceAlbedoImpl->setData(ms);
}

void MU71Atmosphere::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceAlbedo>);
    iIceAlbedoImpl = &Module::getImplementation<IIceAlbedo>();
    tryConfigure(iIceAlbedoImpl);
}

ConfigMap MU71Atmosphere::getConfiguration() const { return {}; }

MU71Atmosphere::HelpMap& MU71Atmosphere::getHelpText(HelpMap& map, bool getAll) { return map; }
MU71Atmosphere::HelpMap& MU71Atmosphere::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    Module::getHelpRecursive<IIceAlbedo>(map, getAll);
    return map;
}

MU71Atmosphere::MU71Atmosphere()
    : iIceAlbedoImpl(nullptr)
    , tsurfAccessor(getStore())
    , hsnowAccessor(getStore())
    , ciceAccessor(getStore())
    , iceAlbedoAccessor(getStore())
    , icePenSWAccessor(getStore())
    , q_sw(monthlyCubicBSpline(swTable))
    , q_lw(monthlyCubicBSpline(lwTable))
    , q_sh(monthlyCubicBSpline(shTable))
    , q_lh(monthlyCubicBSpline(lhTable))
{
}

void MU71Atmosphere::update(const Nextsim::TimestepTime& tst)
{
    iIceAlbedoImpl->update(tst);

    dayOfYear = tst.start.gmtime()->tm_yday;
    isLeap = ((tst.start.gmtime()->tm_year % 4 == 0) && (tst.start.gmtime()->tm_year % 100 != 0))
        || (tst.start.gmtime()->tm_year % 400 == 0);

    // always do update on the host because it currently not worth the effort to port
    // monthlyCubicBSpline
    HField& qow = qowAccessor.getHostRW();
    VField& vwind = vwindAccessor.getHostRW();
    HField& snow = snowAccessor.getHostRW();
    HField& dqia_dt = dqia_dtAccessor.getHostRW();
    HField& evap = evapAccessor.getHostRW();
    HField& rain = rainAccessor.getHostRW();
    HField& subl = sublAccessor.getHostRW();
    HField& qia = qiaAccessor.getHostRW();
    HField& penSW = penSWAccessor.getHostRW();
    UField& uwind = uwindAccessor.getHostRW();
    const AdvectedField& cice = ciceAccessor.getHostRO();
    const AdvectedField& hsnow = hsnowAccessor.getHostRO();
    const AdvectedField& tsurf = tsurfAccessor.getHostRO();
    const HField& iceAlbedo = iceAlbedoAccessor.getHostRO();
    const HField& icePenSW = icePenSWAccessor.getHostRO();

    overElements([&](size_t i) {
        const double Tsurf_K = kelvin(tsurf[i]);

        double sw_in = convFactor * q_sw(dayOfYear, isLeap);
        const double hs = cice[i] > 0 ? hsnow[i] / cice[i] : 0.;
        const double albedoValue = iceAlbedo[i];
        const double i0 = icePenSW[i];
        double qsw = -sw_in * (1. - albedoValue) * (1. - i0);
        penSW[i] = sw_in * (1. - albedoValue) * i0;
        qia[i] = -convFactor
                * (q_sh(dayOfYear, isLeap) + q_lh(dayOfYear, isLeap) + q_lw(dayOfYear, isLeap))
            // LW is tabulated + black body radiation
            + Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 4) + qsw;

        // Just the derivative of the black body radiation
        dqia_dt[i] = 4. * Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 3);

        // Only snowfall if we're not melting
        if ((hs > 0 && tsurf[i] < 0.) || (hs == 0 && tsurf[i] < -Ice::s * Water::mu))
            snow[i] = snowfall();
        else
            snow[i] = 0.;

        // Not needed/specified by M&U '71
        qow[i] = 0.;
        subl[i] = 0.;
        rain[i] = 0.;
        evap[i] = 0.;
        uwind[i] = 0.;
        vwind[i] = 0.;
    });
}

// Snowfall according to M&U '71 (in m/s water equivalent)
double MU71Atmosphere::snowfall()
{
    double const conversionFactor = Ice::rhoSnow / (24. * 3600.);

    // Snowfall rate depends on these dates
    int apr30 = 31 + 28 + 31 + 30;
    if (isLeap)
        ++apr30;

    const int may31 = apr30 + 31;
    const int aug20 = may31 + 30 + 31 + 20;
    const int oct31 = aug20 + 11 + 30 + 31;
    const int dec31 = oct31 + 30 + 31;

    // Snowfall rate in winter: "a linear increase of 5 cm from November 1 to April 30
    const double winterRate = 5e-2 * conversionFactor / double(apr30 + dec31 - oct31);

    if (dayOfYear <= apr30) {
        return winterRate;
    } else if (dayOfYear <= may31) {
        // "an additional 5 cm during the month of May"
        return 5e-2 * conversionFactor / 31.;
    } else if (dayOfYear < aug20) {
        return 0.;
    } else if (dayOfYear <= oct31) {
        /* "a linear accumulation of 30 cm between August 20 and October 30"
         * They don't say anything about October 31 - so I assume they meant 31, not 30 in the quote
         * above. */
        return 30e-2 * conversionFactor / double(oct31 - aug20 + 1);
    } else {
        return winterRate;
    }
}

}
