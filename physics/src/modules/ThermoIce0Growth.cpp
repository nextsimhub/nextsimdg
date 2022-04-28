/*!
 * @file ThermoIce0Growth.cpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0Growth.hpp"

#include "include/IceGrowth.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/ModelArray.hpp"
#include "include/constants.hpp"

namespace Nextsim {

void ThermoIce0Growth::update(const TimestepTime& tsTime)
{
    overElements(std::bind(&ThermoIce0Growth::calculateElement, this, std::placeholders::_1, std::placeholders::_2), tsTime);

}

void ThermoIce0Growth::calculateElement(size_t i, const TimestepTime& tst)
{
    static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
    static const double bulkLHFusionIce = Water::Lf * Ice::rho;

    double remainingFlux = qic[i] - qia[i];
    // Top melt. Melting rate is non-positive.
    double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow;
    snowMelt[i] = snowMeltRate * tst.step;
    double snowSublRate = sublim[i] / Ice::rhoSnow;
    double nowSnow = hsnow[i] + (snowMeltRate - snowSublRate) * tst.step;
    // Use excess flux to melt ice. Non-positive value
    double excessIceMelt = std::min(nowSnow, 0.) * bulkLHFusionSnow / bulkLHFusionIce;
    // With the excess flux noted, clamp the snow thickness to a minimum of zero.
    hsnow[i] = std::max(nowSnow, 0.);
    // Then add snowfall back on top
    hsnow[i] += snowfall[i] * tst.step / Ice::rhoSnow;

    // Bottom melt or growth
    double iceBottomChange = (qic[i] - qio[i]) * tst.step / bulkLHFusionIce;
    // Total thickness change
    deltaHi[i] = excessIceMelt + iceBottomChange;
    hice[i] += deltaHi[i];

    // Amount of melting (only) at the top and bottom of the ice
    topMelt[i] = std::min(excessIceMelt, 0.);
    botMelt[i] = std::min(iceBottomChange, 0.);

    // Snow to ice conversion
    double iceDraught = (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / Water::rhoOcean;

    if (doFlooding && iceDraught > hice[i]) {
        double snowDraught = iceDraught - hice[i];
        snowToIce[i] = snowDraught;
        hsnow[i] -= snowDraught * Ice::rho / Ice::rhoSnow;
        hice[i] = iceDraught;
    } else {
        snowToIce[i] = 0;
    }

    // Melt all ice if it is below minimum threshold
    if (hice[i] < IceGrowth::minimumIceThickness()) {
        if (deltaHi[i] < 0) {
            double scaling = oldHi[i] / deltaHi[i];
            topMelt[i] *= scaling;
            botMelt[i] *= scaling;
        }

        // No snow was converted to ice
        snowToIce[i] = 0.;

        // Change in thickness is all of the old thickness
        deltaHi[i] = -oldHi[i];

        // The ice-ocean flux includes all the latent heat
        qio[i]
            += hice[i] * bulkLHFusionIce / tst.step + hsnow[i] * bulkLHFusionSnow / tst.step;

        // No ice, no snow and the surface temperature is the melting point of ice
        hice[i] = 0.;
        hsnow[i] = 0.;
        tice.zIndexAndLayer(i, 0) = Module::getImplementation<IFreezingPoint>()(sss[i]);
    }
}
} /* namespace Nextsim */
