/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/IceMinima.hpp"
#include "include/KernelAlternatives.hpp"
#include "include/ModelArray.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double ThermoIce0::kappa_s;
static const double k_sDefault = 0.3096;
const double ThermoIce0::freezingPointIce = -Water::mu * Ice::s;

ThermoIce0::ThermoIce0()
    : IIceThermodynamics()
    , snowMeltAccessor(getStore(), RO, ModelArray::Type::H)
    , topMeltAccessor(getStore(), RO, ModelArray::Type::H)
    , botMeltAccessor(getStore(), RO, ModelArray::Type::H)
    , qicAccessor(getStore(), RW, ModelArray::Type::H)
{
    snowMeltAccessor.getHostRW().resize();
    topMeltAccessor.getHostRW().resize();
    botMeltAccessor.getHostRW().resize();
    qicAccessor.getHostRW().resize();
}

void ThermoIce0::update(const TimestepTime& tsTime)
{
    IIceThermodynamics::update(tsTime);

    auto execSpace = DefaultExecutionSpace();
    auto& tsurf = tsurfAccessor.getAutoRW(execSpace);
    auto& deltaHi = deltaHiAccessor.getAutoRW(execSpace);
    auto& hice = hiceAccessor.getAutoRW(execSpace);
    auto& hsnow = hsnowAccessor.getAutoRW(execSpace);
    auto& cice = ciceAccessor.getAutoRW(execSpace);
    auto& qow = qowAccessor.getAutoRW(execSpace);
    auto& qswBase = qswBaseAccessor.getAutoRW(execSpace);
    auto& snowToIce = snowToIceAccessor.getAutoRW(execSpace);
    auto& topMelt = topMeltAccessor.getAutoRW(execSpace);
    auto& snowMelt = snowMeltAccessor.getAutoRW(execSpace);
    auto& botMelt = botMeltAccessor.getAutoRW(execSpace);
    auto& qic = qicAccessor.getAutoRW(execSpace);
    const auto& sublim = sublimAccessor.getAutoRO(execSpace);
    const auto& tf = tfAccessor.getAutoRO(execSpace);
    const auto& snowfall = snowfallAccessor.getAutoRO(execSpace);
    const auto& dQia_dt = dQia_dtAccessor.getAutoRO(execSpace);
    const auto& qio = qioAccessor.getAutoRO(execSpace);
    const auto& penSw = penSwAccessor.getAutoRO(execSpace);
    const auto& qia = qiaAccessor.getAutoRO(execSpace);

    // static members can not be captured directly
    const bool doFlooding = ThermoIce0::doFlooding;
    const double freezingPointIce = ThermoIce0::freezingPointIce;
    const double kappa_s = ThermoIce0::kappa_s;
    const double cMin = IceMinima::c();
    const double hMin = IceMinima::h();
    const double dt = tsTime.step.seconds();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
        static const double bulkLHFusionIce = Water::Lf * Ice::rho;

        // If there is too little ice, do nothing and zero out the computed arrays
        if (hice[i] <= hMin || cice[i] <= cMin) {
            deltaHi[i] = 0.;
            snowToIce[i] = 0.;
            snowMelt[i] = 0.;
            qswBase[i] = 0.;

            tsurf[i] = freezingPointIce;

            // Add to open water flux, since cice will be set to zero
            qow[i] += (hice[i] * bulkLHFusionIce + hsnow[i] * bulkLHFusionSnow) / dt;
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
        tice_i = Utils::min(meltingLimit, tice_i);

        // Top melt. Melting rate is non-positive.
        double snowMeltRate = Utils::min(-remainingFlux, 0.) / bulkLHFusionSnow;
        snowMelt[i] = snowMeltRate * dt;
        double snowSublRate = sublim[i] / Ice::rhoSnow;
        double nowSnow = hs + (snowMeltRate - snowSublRate) * dt;
        // Use excess flux to melt ice. Non-positive value
        double excessIceMelt = Utils::min(nowSnow, 0.) * bulkLHFusionSnow / bulkLHFusionIce;
        // With the excess flux noted, clamp the snow thickness to a minimum of zero.
        hs = Utils::max(nowSnow, 0.);

        // Bottom melt or growth
        double iceBottomChange = (qic[i] - qio[i]) * dt / bulkLHFusionIce;
        // Total thickness change
        deltaHi[i] = excessIceMelt + iceBottomChange;
        hi += deltaHi[i];

        // Then add snowfall back on top if there's still ice
        if (hi > 0.)
            hs += snowfall[i] * dt / Ice::rhoSnow;

        // Amount of melting (only) at the top and bottom of the ice
        topMelt[i] = Utils::min(excessIceMelt, 0.);
        botMelt[i] = Utils::min(iceBottomChange, 0.);
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
        if (hi < hMin) {
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
            qow[i] += cice[i] * (hi * bulkLHFusionIce + hs * bulkLHFusionSnow) / dt;

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
    });
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
                             { "snow_melt", snowMeltAccessor.getHostRO() },
                             { "top_melt", topMeltAccessor.getHostRO() },
                             { "bottom_melt", botMeltAccessor.getHostRO() },
                             { "Q_ic", qicAccessor.getHostRO() },
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

    snowMeltAccessor.getHostRW().resize();
    topMeltAccessor.getHostRW().resize();
    botMeltAccessor.getHostRW().resize();
    qicAccessor.getHostRW().resize();
}

} /* namespace Nextsim */
