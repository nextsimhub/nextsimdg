/*!
 * @file ThermoIce0.cpp
 *
 * @date 29 Jan 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/IceGrowth.hpp"
#include "include/IceMinima.hpp"
#include "include/ModelArray.hpp"
#include "include/NZLevels.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double ThermoIce0::kappa_s;
static const double k_sDefault = 0.3096;
const double ThermoIce0::freezingPointIce = -Water::mu * Ice::s;
const size_t ThermoIce0::nZLevels = 1;

ThermoIce0::ThermoIce0()
    : IIceThermodynamics()
    , botMelt(ModelArray::Type::H)
    , qic(ModelArray::Type::H)
    , topMelt(ModelArray::Type::H)
{
}

void ThermoIce0::update(const TimestepTime& tsTime)
{
    overElements(std::bind(&ThermoIce0::calculateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tsTime);
}

static const std::map<int, std::string> keyMap = {
    { ThermoIce0::KS_KEY, IIceThermodynamics::getKappaSConfigKey() },
};

void ThermoIce0::configure()
{
    kappa_s = Configured::getConfiguration(keyMap.at(KS_KEY), k_sDefault);
    NZLevels::set(nZLevels);
}

ModelState ThermoIce0::getStateRecursive(const OutputSpec& os) const
{
    ModelState state = { {},
        {
            { keyMap.at(KS_KEY), kappa_s },
        } };
    return os ? state : ModelState();
}

ThermoIce0::HelpMap& ThermoIce0::getHelpText(HelpMap& map, bool getAll)
{
    map["ThermoIce0"] = {
        { keyMap.at(KS_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(k_sDefault), "W K⁻¹ m⁻¹", "Thermal conductivity of snow." },
    };
    return map;
}
ThermoIce0::HelpMap& ThermoIce0::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

void ThermoIce0::setData(const ModelState::DataMap& ms)
{
    IIceThermodynamics::setData(ms);

    topMelt.resize();
    botMelt.resize();
    qic.resize();
}

void ThermoIce0::calculateElement(size_t i, const TimestepTime& tst)
{
    // If there is too little ice, do nothing and zero out the computed arrays
    if (hice[i] == 0. || cice[i] == 0.) {
        deltaHi[i] = 0.;
        snowToIce[i] = 0.;

        return;
    }

    static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
    static const double bulkLHFusionIce = Water::Lf * Ice::rho;

    // Semtner's fudge factors for the zero-layer model
    constexpr double beta = 0.4;
    constexpr double gamma = 1.065;

    // Slab thickness and old values
    double hi = hice[i] / cice[i];
    double hs = hsnow[i] / cice[i];
    const double oldHi = hi;

    // Create a reference to the local updated Tice value here to avoid having
    // to write the array access expression out in full every time
    double& tice_i = tice.zIndexAndLayer(i, 0);
    double k_lSlab = kappa_s * Ice::kappa / (kappa_s * hi + Ice::kappa * hs) * gamma;
    qic[i] = k_lSlab * (tf[i] - tice.zIndexAndLayer(i, 0));
    double remainingFlux = qic[i] - (qia[i] + (1. - beta) * penSw[i]);
    tice_i += remainingFlux / (k_lSlab + dQia_dt[i]);

    // Clamp the temperature of the ice to a maximum of the melting point
    // of ice or snow
    double meltingLimit = (hs > 0) ? 0 : freezingPointIce;
    tice_i = std::min(meltingLimit, tice_i);

    // Top melt. Melting rate is non-positive.
    double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow;
    snowMelt[i] = snowMeltRate * tst.step;
    double snowSublRate = sublim[i] / Ice::rhoSnow;
    hs += (snowMeltRate - snowSublRate) * tst.step;
    // Use excess flux to melt ice. Non-positive value
    double excessIceMelt = std::min(hs, 0.) * bulkLHFusionSnow / bulkLHFusionIce;
    // With the excess flux noted, clamp the snow thickness to a minimum of zero.
    hs = std::max(hs, 0.);

    // Bottom melt or growth
    double iceBottomChange = (qic[i] - qio[i]) * tst.step / bulkLHFusionIce;
    // Total thickness change
    deltaHi[i] = excessIceMelt + iceBottomChange;
    hi += deltaHi[i];

    // Then add snowfall back on top if there's still ice
    if (hi > 0.)
        hs += snowfall[i] * tst.step / Ice::rhoSnow;

    // Amount of melting (only) at the top and bottom of the ice
    topMelt[i] = std::min(excessIceMelt, 0.);
    botMelt[i] = std::min(iceBottomChange, 0.);
    // Snow to ice conversion
    double iceDraught = (hi * Ice::rho + hs * Ice::rhoSnow) / Water::rhoOcean;

    if (doFlooding && iceDraught > hi) {
        double snowDraught = iceDraught - hi;
        snowToIce[i] = snowDraught;
        hs -= snowDraught * Ice::rho / Ice::rhoSnow;
        hi = iceDraught;
    } else {
        snowToIce[i] = 0;
    }

    // Melt all ice if it is below minimum threshold
    if (hi < IceMinima::h()) {
        if (deltaHi[i] < 0) {
            double scaling = oldHi / deltaHi[i];
            topMelt[i] *= scaling;
            botMelt[i] *= scaling;
        }

        // No snow was converted to ice
        snowToIce[i] = 0.;

        // Change in thickness is all of the old thickness
        deltaHi[i] = -oldHi;

        // The ice-ocean flux includes all the latent heat
        qio[i] += hi * bulkLHFusionIce / tst.step + hs * bulkLHFusionSnow / tst.step;

        // No ice, no snow and the surface temperature is the melting point of ice
        cice[i] = 0.;
        hi = 0.;
        hs = 0.;
        tice.zIndexAndLayer(i, 0) = celsius(Ice::Tm);
    }

    // Return the cell averaged values
    hice[i] = hi * cice[i];
    hsnow[i] = hs * cice[i];
}

size_t ThermoIce0::getNZLevels() const { return nZLevels; }
} /* namespace Nextsim */
