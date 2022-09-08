//
// Created by Einar Ólason on 01/09/2022.
//

#include "include/MonthlyFluxes.hpp"
#include "include/TableLookup.hpp"

namespace Nextsim {

void MonthlyFluxes::update(const Nextsim::TimestepTime& tst)
{
    overElements(std::bind(&MonthlyFluxes::calculateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

void MonthlyFluxes::calculateElement(size_t i, const TimestepTime& tst)
{
    aoState.update(tst);

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

    IIceAlbedo* iIceAlbedoImpl = &Module::getImplementation<IIceAlbedo>();
    double albedoValue = iIceAlbedoImpl->albedo(tice.zIndexAndLayer(i, 0), h_snow_true[i]);

    // We can set ModelArray to double, which sets all the values.
    qia = (1. - albedoValue) * q_sw + q_lw + q_sh + q_lh;
    qio = 1.5 * convFactor / 12.; // Division by 12 as it's a early average
    qow = 0.;
    subl = 0.;
    // Just the derivative of the black body radiation
    dqia_dt = 4. * Ice::epsilon * PhysicalConstants::sigma * std::pow(Tsurf_K, 3);
}

}