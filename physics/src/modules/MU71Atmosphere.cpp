//
// Created by Einar Ólason on 01/09/2022.
//

#include "include/MU71Atmosphere.hpp"
#include "include/IIceAlbedoModule.hpp"
#include "include/TableLookup.hpp"

namespace Nextsim {

void MU71Atmosphere::configure()
{
    iIceAlbedoImpl = &Module::getImplementation<IIceAlbedo>();
    tryConfigure(iIceAlbedoImpl);
}

void MU71Atmosphere::update(const Nextsim::TimestepTime& tst)
{
    iIceAlbedoImpl->setTime(tst.start);
    overElements(std::bind(&MU71Atmosphere::calculateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

void MU71Atmosphere::calculateElement(size_t i, const TimestepTime& tst)
{
    // Monthly fluxes from Maykut and Untersteiner (1971)
    const std::vector<double> swTable
        = { 0.00, 0.00, 1.90, 9.99, 17.7, 19.2, 13.6, 9.00, 3.70, 0.40, 0.00, 0.00 };
    const std::vector<double> lwTable
        = { 10.4, 10.3, 10.3, 11.6, 15.1, 18.0, 19.1, 18.7, 16.5, 13.9, 11.2, 10.9 };
    const std::vector<double> shTable
        = { 1.18, 0.76, 0.72, 0.29, -.45, -.39, -.30, -.40, -.17, 0.10, 0.56, 0.79 };
    const std::vector<double> lhTable
        = { 0.00, -.02, -.03, -.09, -.46, -.70, -.64, -.66, -.39, -.19, -.01, -.01 };
    //       Jan,  Feb,  Mar,  Apr,  Mai,  Jun,  Jul,  Aug,  Sept, Oct,  Nov,  Dec

    // Conversion factor from kcal/cm^2/month to W/m^2
    const double convFactor = 4.184e7 / (365.2425 / 12. * 24. * 3600.);

    const double dayOfYear = tst.start.gmtime()->tm_yday;
    const bool isLeap
        = ((tst.start.gmtime()->tm_year % 4 == 0) && (tst.start.gmtime()->tm_year % 100 != 0))
        || (tst.start.gmtime()->tm_year % 400 == 0);

    // Just tabulated values
    const double q_sw = -convFactor * TableLookup::monthlyLinearLUT(swTable, dayOfYear, isLeap);
    const double q_sh = convFactor * TableLookup::monthlyLinearLUT(shTable, dayOfYear, isLeap);
    const double q_lh = convFactor * TableLookup::monthlyLinearLUT(lhTable, dayOfYear, isLeap);

    // LW is tabulated + black body radiation
    const double Tsurf_K = tice.zIndexAndLayer(i, 0) + PhysicalConstants::Tt;
    const double q_lw = -convFactor * TableLookup::monthlyLinearLUT(lwTable, dayOfYear, isLeap)
        + Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 4);

    double albedoValue = iIceAlbedoImpl->albedo(tice.zIndexAndLayer(i, 0), h_snow_true[i]);

    // We can set ModelArray to double, which sets all the values.
    qia = (1. - albedoValue) * q_sw + q_lw + q_sh + q_lh;
    // Just the derivative of the black body radiation
    dqia_dt = 4. * Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 3);

    // Only snowfall if we're not melting
    if (tice.zIndexAndLayer(i, 0) <= 0.)
        snow = snowfall(dayOfYear, isLeap);
    else
        snow = 0.;

    // Not needed/specified by M&U '71
    qow = 0.;
    subl = 0.;
    rain = 0.;
    evap = 0.;
    uwind = 0.;
    vwind = 0.;
}

// Snowfall according to M&U '71 (in m/s water equivalent)
double MU71Atmosphere::snowfall(const double dayOfYear, bool isLeap = false)
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
        return 5e-2 * conversionFactor / 30.;
    } else if (dayOfYear <= aug20) {
        return 0.;
    } else if (dayOfYear <= oct31) {
        /* "a linear accumulation of 30 cm between August 20 and October 30"
         * They don't say anything about October 31 - so I assume they meant 31, not 30 in the quote
         * above. */
        return 30e-2 * conversionFactor / double(oct31 - aug20);
    } else {
        return winterRate;
    }
}

}