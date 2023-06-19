/*!
 * @file ThermoWinton.cpp
 *
 * @date Sep 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoWinton.hpp"
#include "include/IceMinima.hpp"
#include "include/NZLevels.hpp"

#include "include/constants.hpp"

#include <cmath>

namespace Nextsim {

const size_t ThermoWinton::nLevels = 3;
double ThermoWinton::kappa_s;
double ThermoWinton::i0;
static const double k_sDefault = 0.3096;
static const double i0_default = 0.17;

const double ThermoWinton::cVol = Ice::cp * Ice::rho; // bulk heat capacity of ice
const double ThermoWinton::seaIceTf = -Water::mu * Ice::s;
bool ThermoWinton::doFlooding = true;

ThermoWinton::ThermoWinton()
    : IIceThermodynamics()
    , snowMelt(ModelArray::Type::H)
    , topMelt(ModelArray::Type::H)
    , botMelt(ModelArray::Type::H)
    , oldHi(getProtectedArray())
    , sw_in(getProtectedArray())
    , subl(getSharedArray())
{
}

template <>
const std::map<int, std::string> Configured<ThermoWinton>::keyMap = {
    { ThermoWinton::KS_KEY, IIceThermodynamics::getKappaSConfigKey() },
    { ThermoWinton::I0_KEY, "nextsim_thermo.I_0" },
    { ThermoWinton::FLOODING_KEY, "nextsim_thermo.doFlooding" },
};

void ThermoWinton::configure()
{
    kappa_s = Configured::getConfiguration(keyMap.at(KS_KEY), k_sDefault);
    i0 = Configured::getConfiguration(keyMap.at(I0_KEY), i0_default);
    doFlooding = Configured::getConfiguration(keyMap.at(FLOODING_KEY), doFlooding);
    NZLevels::set(nLevels);
}

ModelState ThermoWinton::getStateRecursive(const OutputSpec& os) const
{
    ModelState state = { {},
        {
            { keyMap.at(KS_KEY), kappa_s },
            { keyMap.at(I0_KEY), i0 },
        } };
    return os ? state : ModelState();
}

ThermoWinton::HelpMap& ThermoWinton::getHelpText(HelpMap& map, bool getAll)
{
    map["ThermoWinton"] = {
        { keyMap.at(KS_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(k_sDefault),
            "W K⁻¹ m⁻¹", "Thermal conductivity of snow." },
        { keyMap.at(I0_KEY), ConfigType::NUMERIC, { "0", "1" }, std::to_string(i0_default),
            "unitless", "Optical albedo of liquid water." },
    };
    return map;
}
ThermoWinton::HelpMap& ThermoWinton::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

void ThermoWinton::setData(const ModelState::DataMap& state)
{
    IIceThermodynamics::setData(state);

    snowMelt.resize();
    topMelt.resize();
    botMelt.resize();
    snowToIce.resize();

    // The Winton scheme requires three temperature levels in the ice
    if (tice0.data().size() != nLevels * hice.data().size()) {
        double actualLevels = static_cast<double>(tice0.data().size()) / hice.data().size();
        throw std::length_error(std::string("The inferred number of ice temperature levels is ")
            + std::to_string(actualLevels) + " when the Winton ice thermodynamics scheme expects "
            + std::to_string(nLevels));
    }
}

void ThermoWinton::update(const TimestepTime& tst)
{
    overElements(std::bind(&ThermoWinton::calculateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

size_t ThermoWinton::getNZLevels() const { return nLevels; }

void ThermoWinton::calculateElement(size_t i, const TimestepTime& tst)
{

    // Don't do anything if there is no ice
    if (cice[i] <= 0 || hice[i] <= 0) {

        snowToIce[i] = 0;

        deltaHi[i] = 0;
        hice[i] = 0;
        hsnow[i] = 0;

        tice.zIndexAndLayer(i, 0) = seaIceTf;
        tice.zIndexAndLayer(i, 1) = seaIceTf;
        tice.zIndexAndLayer(i, 2) = seaIceTf;

        return;
    }

    static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
    static const double bulkLHFusionIce = Water::Lf * Ice::rho;

    double tSurf = tice0.zIndexAndLayer(i, 0); // surface temperature
    double tUppr = tice0.zIndexAndLayer(i, 1); // upper layer temperature
    double tLowr = tice0.zIndexAndLayer(i, 2); // lower layer temperature
    double tBott = tf[i]; // freezing point of (local) seawater

    double dt = tst.step.seconds();

    double surfMelt = 0; // surface melting mass loss
    // Calculate temperatures by solving the heat conduction equation
    calculateTemps(tSurf, tUppr, tLowr, surfMelt, i, dt);

    // The ratio of ΔΗ_f T_f / c_p,ice is used a lot. Units are K²
    const static double dHfTf_cp = Water::Lf * seaIceTf / Ice::cp;

    // Thickness changes
    // ice
    double h1 = hice[i] / 2;
    double h2 = hice[i] / 2;
    // Eqs. (1) and (25) - but I 've multiplied them with \rho_i (hence cVol), because it's missing
    // in the paper
    double e1 = cVol * (tUppr - seaIceTf) - bulkLHFusionIce * (1 - seaIceTf / tUppr);
    double e2 = cVol * (tLowr - seaIceTf) - bulkLHFusionIce;

    double& hs = hsnow[i];
    // snow
    hs += snowfall[i] / Ice::rhoSnow * dt;
    //    double accumulatedSnowThickness = snowfall[i] / Ice::rhoSnow * dt;

    // sublimation
    // 4 cases
    const double& subli = subl[i];
    double deltaSnow = subli * dt / Ice::rhoSnow;
    double deltaIce1 = (subli * dt - hs * Ice::rhoSnow) / Ice::rho;
    double deltaIce2 = deltaIce1 - h1;
    if (deltaSnow <= hs) {
        // sublimation is less than or equal to the mass of snow
        hs -= deltaSnow;
    } else if (deltaIce1 <= h1) {
        // sublimation minus sublimed snow is less than or equal to half the
        // ice thickness
        h1 -= deltaIce1;
        hs = 0;
    } else if (deltaIce2 <= h2) {
        // sublimation minus sublimed snow is greater than half the ice
        // thickness, but not all of it
        h2 -= deltaIce2;
        h1 = 0;
        hs = 0;
    } else {
        // the snow and ice sublimates
        double oceanEvapError = (deltaIce2 - h2) * Ice::rho / Water::rhoOcean;
        // TODO: log the error
        h2 = 0;
        h1 = 0;
        hs = 0;
    }
    // Sublimated ice counts as top melt
    topMelt[i] = std::max(0., h1 + h2 - hice[i]); // (23)

    // Bottom melt/freezing
    double meltBottom = (qio[i] - 4 * Ice::kappa * (tBott - tLowr) / hice[i]) * dt;
    snowMelt[i] = 0;
    if (meltBottom <= 0.) {
        // Freezing
        // Eq. (25) - but I 've multiplied them with \rho_i (hence cVol), because it's missing in
        // the paper
        double eBot = cVol * (tBott - seaIceTf) - bulkLHFusionIce;
        deltaIce2 = meltBottom / eBot; // (24)
        tLowr = (deltaIce2 * tBott + h2 * tLowr) / (deltaIce2 + h2); // (26)
        h2 += deltaIce2;
    } else {
        // Melting
        // Eqs. (31)-(32) with added division with \rho_i (and \rho_s for 32)
        deltaIce2 = -std::min(-meltBottom / e2, h2);
        deltaIce1 = -std::min(std::max(-(meltBottom + e2 * h2) / e1, 0.), h1);
        snowMelt[i] = -std::min(
            std::max((meltBottom + e2 * h2 + e1 * h1) / bulkLHFusionSnow, 0.), hsnow[i]);

        // If everything melts we need to put heat back into the ocean
        if (h2 + h1 + hs - deltaIce2 - deltaIce1 - snowMelt[i] <= 0.) {
            // (34) - with added multiplication of rhoi and rhos and division with dt
            qio[i] -= std::max(meltBottom - bulkLHFusionSnow * hs + e1 * h1 + e2 * h2, 0.) / dt;
        }

        hs += snowMelt[i];
        h1 += deltaIce1;
        h2 += deltaIce2;
        botMelt[i] += deltaIce1 + deltaIce2;
    }

    // Melting at the surface
    // Do we really need an assertion here?
    // assert(surfMelt >= 0);
    // Eqs. (27)-(29) with division of \rho_i and \rho_s
    snowMelt[i] -= std::min(surfMelt * dt / bulkLHFusionSnow, hs);
    deltaIce1 = -std::min(std::max(-(surfMelt * dt - bulkLHFusionSnow * hs) / e1, 0.), h1);
    deltaIce2
        = -std::min(std::max(-(surfMelt * dt - bulkLHFusionSnow * hs + e1 * h1) / e2, 0.), h2);

    // If everything melts we need to put heat back into the ocean
    // Eq (30) - with multiplication of rhoi and rhos and division with dt
    if (h2 + h1 + hs - deltaIce2 - deltaIce1 - snowMelt[i] <= 0.) {
        qio[i] -= std::max(surfMelt * dt - bulkLHFusionSnow * hs + e1 * h1 + e2 * h2, 0.) / dt;
    }

    hs += snowMelt[i];
    h1 += deltaIce1;
    h2 += deltaIce2;
    topMelt[i] += deltaIce1 + deltaIce2;

    // Snow to ice conversion
    double freeboard
        = (hice[i] * (Water::rhoOcean - Ice::rho) - hs * Ice::rhoSnow) / Water::rhoOcean;
    if (doFlooding && freeboard < 0.) {
        hs += std::min(freeboard * Ice::rho / Ice::rhoSnow, 0.); // (35) using +=
        deltaIce1 = std::max(-freeboard, 0.); // (36)
        double f1 = 1 - deltaIce1 / (deltaIce1 + h1); // Fraction of new ice in the upper layer
        double tBar = f1 * (tUppr + dHfTf_cp / tUppr) + (1 - f1) * seaIceTf; // (39)
        tUppr = (tBar - std::sqrt(tBar * tBar - 4 * dHfTf_cp)) / 2; // (38)
        h1 += deltaIce1;
        snowToIce[i] += deltaIce1;
    }

    // Add up the half-layer thicknesses
    double& hi = hice[i];
    hi = h1 + h2;
    // Adjust the temperatures to evenly divide the ice
    if (h2 > h1) {
        // Lower layer ice is added to the upper layer
        double f1 = h1 / hi * 2;
        double tBar = f1 * (tUppr + dHfTf_cp / tUppr) + (1 - f1) * tLowr; // (39)
        // The upper layer temperature changes
        tUppr = (tBar - std::sqrt(tBar * tBar - 4 * dHfTf_cp)) / 2; // (38)
    } else {
        // Upper layer ice is added to the lower layer
        double f1 = (2 * h1 - hi) / hi;
        // Lower layer temperature changes
        tLowr = f1 * (tUppr + dHfTf_cp / tUppr) + (1 - f1) * tLowr; // (40)
        // Melt from top and bottom if the lower layer temperature is too high
        if (tLowr > seaIceTf) {
            double deltaMelt = hi / 4 * Ice::cp * (tLowr - seaIceTf) * tUppr
                / (Ice::Lf * tUppr + (Ice::cp * tUppr - Ice::Lf) * (seaIceTf - tUppr));
            topMelt[i] -= deltaMelt;
            botMelt[i] -= deltaMelt;
            hi -= 2 * deltaMelt;
            tLowr = seaIceTf;
        }
    }
    deltaHi[i] = hi - oldHi[i];

    // Remove very small ice thickness
    if (hi < IceMinima::h()) {
        // (30) - with multiplication of rhoi and rhos and division with dt
        qio[i] -= (-bulkLHFusionSnow * hs + (e1 + e2) * hi / 2) / dt;

        if (deltaHi[i] < 0) {
            topMelt[i] *= oldHi[i] / deltaHi[i];
            botMelt[i] *= oldHi[i] / deltaHi[i];
        }
        snowToIce[i] = 0;

        deltaHi[i] = -oldHi[i];
        hi = 0;
        hs = 0;
        tSurf = seaIceTf;
        tUppr = seaIceTf;
        tLowr = seaIceTf;
    }

    tice.zIndexAndLayer(i, 0) = tSurf;
    tice.zIndexAndLayer(i, 1) = tUppr;
    tice.zIndexAndLayer(i, 2) = tLowr;
}

void ThermoWinton::calculateTemps(
    double& tSurf, double& tUppr, double& tLowr, double& mSurf, size_t i, double dt)
{
    /*
     * Internal temperatures
     * Following Winton (2000) we start by calculating the temperature - first T1, then Tsurf, and
     * finally T2 Numers in parentheses refer to equations in the paper
     */

    double& hi = hice[i];
    double tBase = tf[i]; // Freezing point of seawater with the local salinity
    double tMelt = (hsnow[i] > 0) ? 0 : seaIceTf; // Melting point at the surface

    // First some coefficients based on temperatures from the previous time step
    double k12 = 4 * Ice::kappa * kappa_s / (kappa_s * hi + 4 * Ice::kappa * hsnow[i]); // (5)
    double a = qia[i] - tSurf * dQia_dt[i]; // (7)
    double b = dQia_dt[i]; // (8)
    double k32 = 2 * Ice::kappa / hi; // (10)

    double a1 = hi * cVol / (2 * dt) + k32 * (4 * dt * k32 + hi * cVol) / (6 * dt * k32 + hi * cVol)
        + b * k12 / (k12 + b); // (16)
    double b1 = -hi * (cVol * tUppr + Ice::Lf * Ice::rho * seaIceTf / tUppr) / (2 * dt) - penSw[i]
        - k32 * (4 * dt * k32 * tBase + hi * cVol * tLowr) / (6 * dt * k32 + hi * cVol)
        + a * k12 / (k12 + b); // (17)
    double c1 = hi * Ice::Lf * Ice::rho * seaIceTf / (2 * dt); // (18)

    // Updated surface and mid-ice temperatures
    tUppr = -(b1 + std::sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1); // (21)
    tSurf = (k12 * tUppr - a) / (k12 + b); // (6)
    // Is the surface melting?
    if (tSurf > tMelt) {
        // Recalculate T1 and Tsurf if so
        tSurf = tMelt;
        // apply the change to the *1 parameters, rather than recalculating in full
        a1 += k12 - k12 * b / (k12 + b); // (19) simplified using +=
        b1 -= k12 * tSurf + a * k12 / (k12 + b); // (20) simplified using -=
        tUppr = -(b1 + std::sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1); // (21)

        // Surface melt
        mSurf = k12 * (tUppr - tSurf) - (a + b * tSurf); // (22)
    }

    // update the lower T based on the new value of the upper T (eq 15)
    tLowr = (2 * dt * k32 * (tUppr + 2 * tf[i]) + hi * cVol * tLowr) / (6 * dt * k32 + hi * cVol);
}

} /* namespace Nextsim */
