//
// Created by Einar Ólason on 01/09/2022.
//

#include "include/MU71Atmosphere.hpp"
#include "include/TableLookup.hpp"

namespace Nextsim {

void MU71Atmosphere::update(const Nextsim::TimestepTime& tst)
{
    dayOfYear = tst.start.gmtime()->tm_yday;
    isLeap = ((tst.start.gmtime()->tm_year % 4 == 0) && (tst.start.gmtime()->tm_year % 100 != 0))
        || (tst.start.gmtime()->tm_year % 400 == 0);

    overElements(std::bind(&MU71Atmosphere::calculateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

void MU71Atmosphere::calculateElement(size_t i, const TimestepTime& tst)
{
    // Just tabulated values
    const double q_sw = -convFactor * TableLookup::monthlyLinearLUT(swTable, dayOfYear, isLeap);
    const double q_sh = -convFactor * TableLookup::monthlyLinearLUT(shTable, dayOfYear, isLeap);
    const double q_lh = -convFactor * TableLookup::monthlyLinearLUT(lhTable, dayOfYear, isLeap);

    // LW is tabulated + black body radiation
    const double Tsurf_K = tice.zIndexAndLayer(i, 0) + PhysicalConstants::Tt;
    const double q_lw = -convFactor * TableLookup::monthlyLinearLUT(lwTable, dayOfYear, isLeap)
        + Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 4);

    // We can set ModelArray to double, which sets all the values.
    qsw = q_sw;
    qia = q_lw + q_sh + q_lh;
    // Just the derivative of the black body radiation
    dqia_dt = 4. * Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 3);

    // Only snowfall if we're not melting
    if (tice.zIndexAndLayer(i, 0) <= 0.)
        snow = snowfall();
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