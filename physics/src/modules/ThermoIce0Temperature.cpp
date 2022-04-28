/*!
 * @file ThermoIce0Temperature.cpp
 *
 * @date Apr 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0Temperature.hpp"

namespace Nextsim {
double ThermoIce0Temperature::k_s;
const double ThermoIce0Temperature::freezingPointIce = -Water::mu * Ice::s;

template <>
const std::map<int, std::string> Configured<ThermoIce0Temperature>::keyMap = {
    { ThermoIce0Temperature::KS_KEY, "thermoice0.ks" },
};

void ThermoIce0Temperature::configure()
{
    k_s = Configured::getConfiguration(keyMap.at(KS_KEY), 0.3096);
}

void ThermoIce0Temperature::calculateElement(size_t i, const TimestepTime& tst)
{
    double& tice_i = tice.zIndexAndLayer(i, 0);
    double k_lSlab = k_s * Ice::kappa / (k_s * hice[i] + Ice::kappa * hsnow[i]);
    qic[i] = k_lSlab * (tf[i] - tice0.zIndexAndLayer(i, 0));
    double remainingFlux = qic[i] - qia[i];
    tice_i = tice0.zIndexAndLayer(i, 0) + remainingFlux / (k_lSlab + dqia_dt[i]);

    // Clamp the temperature of the ice to a maximum of the melting point
    // of ice or snow
    double meltingLimit = (hsnow[i] > 0) ? 0 : freezingPointIce;
    tice_i = std::min(meltingLimit, tice_i);
}

} /* namespace Nextsim */
