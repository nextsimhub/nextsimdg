/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/IceGrowth.hpp"
#include "include/IceMinima.hpp"
#include "include/ModelArray.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double ThermoIce0::kappa_s;
static const double k_sDefault = 0.3096;
const double ThermoIce0::freezingPointIce = -Water::mu * Ice::s;

ThermoIce0::ThermoIce0()
    : IIceThermodynamics()
    , snowMelt(ModelArray::Type::H)
    , topMelt(ModelArray::Type::H)
    , botMelt(ModelArray::Type::H)
    , qic(ModelArray::Type::H)
{
}

void ThermoIce0::update(const TimestepTime& tsTime)
{
    IIceThermodynamics::update(tsTime);

    AdvectedField& tsurf = tsurfAccessor.getHostRW();
    HField& deltaHi = deltaHiAccessor.getHostRW();
    AdvectedField& hice = hiceAccessor.getHostRW();
    AdvectedField& hsnow = hsnowAccessor.getHostRW();
    AdvectedField& cice = ciceAccessor.getHostRW();
    HField& qow = qowAccessor.getHostRW();
    HField& qswBase = qswBaseAccessor.getHostRW();
    HField& snowToIce = snowToIceAccessor.getHostRW();
    const HField& sublim = sublimAccessor.getHostRO();
    const HField& tf = tfAccessor.getHostRO();
    const HField& snowfall = snowfallAccessor.getHostRO();
    const HField& dQia_dt = dQia_dtAccessor.getHostRO();
    const HField& qio = qioAccessor.getHostRO();
    const HField& penSw = penSwAccessor.getHostRO();
    const HField& qia = qiaAccessor.getHostRO();

    overElements(
        [&](size_t i, const TimestepTime& tst) {
            static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
            static const double bulkLHFusionIce = Water::Lf * Ice::rho;

            // If there is too little ice, do nothing and zero out the computed arrays
            if (hice[i] <= IceMinima::h() || cice[i] <= IceMinima::c()) {
                deltaHi[i] = 0.;
                snowToIce[i] = 0.;
                snowMelt[i] = 0.;
                qswBase[i] = 0.;

                tsurf[i] = freezingPointIce;

                // Add to open water flux, since cice will be set to zero
                qow[i] += (hice[i] * bulkLHFusionIce + hsnow[i] * bulkLHFusionSnow) / tst.step;
                cice[i] = 0.;
                hice[i] = 0.;
                hsnow[i] = 0.;

                return;
            }

            // Semtner's fudge factors for the zero-layer model
            constexpr double beta = 0.4;
            constexpr double gamma = 1.065;

            double hi = hice[i] / cice[i];
            const double oldHi = hi;
            double hs = hsnow[i] / cice[i];
            // Create a reference to the local updated Tice value here to avoid having
            // to write the array access expression out in full every time
            double& tice_i = tsurf[i];
            double tice0 = tsurf[i];
            double k_lSlab = kappa_s * Ice::kappa / (kappa_s * hi + Ice::kappa * hs) * gamma;
            qic[i] = k_lSlab * (tf[i] - tice0);
            double remainingFlux = qic[i] - (qia[i] + (1. - beta) * penSw[i]);
            tice_i = tice0 + remainingFlux / (k_lSlab + dQia_dt[i]);

            // Clamp the temperature of the ice to a maximum of the melting point
            // of ice or snow
            double meltingLimit = (hs > 0) ? 0 : freezingPointIce;
            tice_i = std::min(meltingLimit, tice_i);

            // Top melt. Melting rate is non-positive.
            double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow;
            snowMelt[i] = snowMeltRate * tst.step;
            double snowSublRate = sublim[i] / Ice::rhoSnow;
            double nowSnow = hs + (snowMeltRate - snowSublRate) * tst.step;
            // Use excess flux to melt ice. Non-positive value
            double excessIceMelt = std::min(nowSnow, 0.) * bulkLHFusionSnow / bulkLHFusionIce;
            // With the excess flux noted, clamp the snow thickness to a minimum of zero.
            hs = std::max(nowSnow, 0.);

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

                // Add the melt flux to open water flux, since cice will be set to zero
                qow[i] += cice[i] * (hi * bulkLHFusionIce + hs * bulkLHFusionSnow) / tst.step;

                // No ice, no snow and the surface temperature is the melting point of ice
                cice[i] = 0.;
                hice[i] = 0.;
                hsnow[i] = 0.;
                tsurf[i] = celsius(Ice::Tm);
            } else {
                // If there is still ice, we need to update the DG ice and snow thicknesses with the
                // new values
                hice[i] = hi * cice[i];
                hsnow[i] = hs * cice[i];
            }
        },
        tsTime);
}

static const std::map<int, std::string> keyMap = {
    { ThermoIce0::KS_KEY, IIceThermodynamics::getKappaSConfigKey() },
};

void ThermoIce0::configure()
{
    kappa_s = Configured::getConfiguration(keyMap.at(KS_KEY), k_sDefault);
}

ConfigMap ThermoIce0::getConfiguration() const { return { { keyMap.at(KS_KEY), kappa_s } }; }

ModelState ThermoIce0::getStateDiagnostic() const
{
    ModelState state = { {
                             { "snow_melt", snowMelt },
                             { "top_melt", topMelt },
                             { "bottom_melt", botMelt },
                             { "Q_ic", qic },
                         },
        getConfiguration() };

    return state.merge(IIceThermodynamics::getStateDiagnostic());
}

ModelState ThermoIce0::getStatePrognostic() const
{
    ModelState state = IIceThermodynamics::getStatePrognostic();
    return state.merge(getConfiguration());
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

    snowMelt.resize();
    topMelt.resize();
    botMelt.resize();
    qic.resize();
}

// void ThermoIce0::calculateElement(size_t i, const TimestepTime& tst)
// {
//     static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
//     static const double bulkLHFusionIce = Water::Lf * Ice::rho;
//
//     // If there is too little ice, do nothing and zero out the computed arrays
//     if (hice[i] <= IceMinima::h() || cice[i] <= IceMinima::c()) {
//         deltaHi[i] = 0.;
//         snowToIce[i] = 0.;
//         snowMelt[i] = 0.;
//         qswBase[i] = 0.;
//
//         tsurf[i] = freezingPointIce;
//
//         // Add to open water flux, since cice will be set to zero
//         qow[i] += (hice[i] * bulkLHFusionIce + hsnow[i] * bulkLHFusionSnow) / tst.step;
//         cice[i] = 0.;
//         hice[i] = 0.;
//         hsnow[i] = 0.;
//
//         return;
//     }
//
//     // Semtner's fudge factors for the zero-layer model
//     constexpr double beta = 0.4;
//     constexpr double gamma = 1.065;
//
//     double hi = hice[i] / cice[i];
//     const double oldHi = hi;
//     double hs = hsnow[i] / cice[i];
//     // Create a reference to the local updated Tice value here to avoid having
//     // to write the array access expression out in full every time
//     double& tice_i = tsurf[i];
//     double tice0 = tsurf[i];
//     double k_lSlab = kappa_s * Ice::kappa / (kappa_s * hi + Ice::kappa * hs) * gamma;
//     qic[i] = k_lSlab * (tf[i] - tice0);
//     double remainingFlux = qic[i] - (qia[i] + (1. - beta) * penSw[i]);
//     tice_i = tice0 + remainingFlux / (k_lSlab + dQia_dt[i]);
//
//     // Clamp the temperature of the ice to a maximum of the melting point
//     // of ice or snow
//     double meltingLimit = (hs > 0) ? 0 : freezingPointIce;
//     tice_i = std::min(meltingLimit, tice_i);
//
//     // Top melt. Melting rate is non-positive.
//     double snowMeltRate = std::min(-remainingFlux, 0.) / bulkLHFusionSnow;
//     snowMelt[i] = snowMeltRate * tst.step;
//     double snowSublRate = sublim[i] / Ice::rhoSnow;
//     double nowSnow = hs + (snowMeltRate - snowSublRate) * tst.step;
//     // Use excess flux to melt ice. Non-positive value
//     double excessIceMelt = std::min(nowSnow, 0.) * bulkLHFusionSnow / bulkLHFusionIce;
//     // With the excess flux noted, clamp the snow thickness to a minimum of zero.
//     hs = std::max(nowSnow, 0.);
//
//     // Bottom melt or growth
//     double iceBottomChange = (qic[i] - qio[i]) * tst.step / bulkLHFusionIce;
//     // Total thickness change
//     deltaHi[i] = excessIceMelt + iceBottomChange;
//     hi += deltaHi[i];
//
//     // Then add snowfall back on top if there's still ice
//     if (hi > 0.)
//         hs += snowfall[i] * tst.step / Ice::rhoSnow;
//
//     // Amount of melting (only) at the top and bottom of the ice
//     topMelt[i] = std::min(excessIceMelt, 0.);
//     botMelt[i] = std::min(iceBottomChange, 0.);
//     // Snow to ice conversion
//     double iceDraught = (hi * Ice::rho + hs * Ice::rhoSnow) / Water::rhoOcean;
//
//     if (doFlooding && iceDraught > hi) {
//         double snowDraught = iceDraught - hi;
//         snowToIce[i] = snowDraught;
//         hs -= snowDraught * Ice::rho / Ice::rhoSnow;
//         hi = iceDraught;
//     } else {
//         snowToIce[i] = 0;
//     }
//
//     // Melt all ice if it is below minimum threshold
//     if (hi < IceMinima::h()) {
//         if (deltaHi[i] < 0) {
//             double scaling = oldHi / deltaHi[i];
//             topMelt[i] *= scaling;
//             botMelt[i] *= scaling;
//         }
//
//         // No snow was converted to ice
//         snowToIce[i] = 0.;
//
//         // Change in thickness is all of the old thickness
//         deltaHi[i] = -oldHi;
//
//         // Add the melt flux to open water flux, since cice will be set to zero
//         qow[i] += cice[i] * (hi * bulkLHFusionIce + hs * bulkLHFusionSnow) / tst.step;
//
//         // No ice, no snow and the surface temperature is the melting point of ice
//         cice[i] = 0.;
//         hice[i] = 0.;
//         hsnow[i] = 0.;
//         tsurf[i] = celsius(Ice::Tm);
//     } else {
//         // If there is still ice, we need to update the DG ice and snow thicknesses with the new
//         // values
//         hice[i] = hi * cice[i];
//         hsnow[i] = hs * cice[i];
//     }
// }

} /* namespace Nextsim */
