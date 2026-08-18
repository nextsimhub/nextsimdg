/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/ThermoWinton.hpp"
#include "include/IceMinima.hpp"

#include "include/KernelAlternatives.hpp"
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif
#include "include/constants.hpp"
#include "include/gridNames.hpp"

#include <cmath>

namespace Nextsim {

FloatType ThermoWinton::kappa_s;
static const FloatType k_sDefault = 0.3096;

const FloatType ThermoWinton::cVol = Ice::cp * Ice::rho; // bulk heat capacity of ice
const FloatType ThermoWinton::seaIceTf = -Water::mu * Ice::s;
bool ThermoWinton::doFlooding = true;

// Names of the additional ice temperature fields
const std::string ThermoWinton::tInteriorName = "tinterior";
const std::string ThermoWinton::tBottomName = "tbottom";

ThermoWinton::ThermoWinton()
    : IIceThermodynamics()
    , tInternalAccessor(getStore(), RO, ModelArray::AdvectionType)
    , tBottomAccessor(getStore(), RO, ModelArray::AdvectionType)
    , snowMeltAccessor(getStore(), RO, ModelArray::Type::H)
    , topMeltAccessor(getStore(), RO, ModelArray::Type::H)
    , botMeltAccessor(getStore(), RO, ModelArray::Type::H)
    , sw_inAccessor(getStore())
    , sublAccessor(getStore())
{
    tInternalAccessor.getHostRW().reinitialize();
    tBottomAccessor.getHostRW().reinitialize();
    snowMeltAccessor.getHostRW().reinitialize();
    topMeltAccessor.getHostRW().reinitialize();
    botMeltAccessor.getHostRW().reinitialize();
    snowToIceAccessor.getHostRW().reinitialize();

#ifdef USE_XIOS
    // Set XIOS field types
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.setPrognosticFieldType(tInteriorName, ModelArray::AdvectionType);
    xiosHandler.setPrognosticFieldType(tBottomName, ModelArray::AdvectionType);
#endif
}

static const std::map<int, std::string> keyMap = {
    { ThermoWinton::KS_KEY, IIceThermodynamics::getKappaSConfigKey() },
    { ThermoWinton::FLOODING_KEY, "nextsim_thermo.doFlooding" },
};

void ThermoWinton::configure()
{
    kappa_s = Configured::getConfiguration(keyMap.at(KS_KEY), k_sDefault);
    doFlooding = Configured::getConfiguration(keyMap.at(FLOODING_KEY), doFlooding);
}

ConfigMap ThermoWinton::getConfiguration() const
{
    return {
        { keyMap.at(KS_KEY), kappa_s },
        { keyMap.at(FLOODING_KEY), doFlooding },
    };
}

ModelState ThermoWinton::getStateDiagnostic() const
{
    ModelState state = { {
                             { "snow_melt", snowMeltAccessor.getHostRO() },
                             { "top_melt", topMeltAccessor.getHostRO() },
                             { "bottom_melt", botMeltAccessor.getHostRO() },
                         },
        getConfiguration() };

    state.merge(getStatePrognostic());
    return state.merge(IIceThermodynamics::getStateDiagnostic());
}

ModelState ThermoWinton::getStatePrognostic() const
{
    ModelState state = { {
                             { tInteriorName, tInternalAccessor.getHostRO() },
                             { tBottomName, tBottomAccessor.getHostRO() },
                         },
        getConfiguration() };

    return state.merge(IIceThermodynamics::getStatePrognostic());
}

ThermoWinton::HelpMap& ThermoWinton::getHelpText(HelpMap& map, bool getAll)
{
    map["ThermoWinton"] = {
        { keyMap.at(KS_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(k_sDefault), "W K⁻¹ m⁻¹", "Thermal conductivity of snow." },
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

    AdvectedField& tInternal = tInternalAccessor.getHostRW();
    AdvectedField& tBottom = tBottomAccessor.getHostRW();
    tInternalAccessor.getHostRW().reinitialize();
    tBottomAccessor.getHostRW().reinitialize();
    snowMeltAccessor.getHostRW().reinitialize();
    topMeltAccessor.getHostRW().reinitialize();
    botMeltAccessor.getHostRW().reinitialize();
    snowToIceAccessor.getHostRW().reinitialize();

    /* If the internal temperature is not in the restart file, then we simply set it to the freezing
     * point of seawater. It's a safe approximation, and it seems the user doesn't really care! */
    try {
        tInternal = state.at(tInteriorName);
    } catch (const std::out_of_range& e) {
        Logged::info("No " + tInteriorName
            + " field in restart file. Setting it to the melting point of ice.\n");
        tInternal = seaIceTf;
    }

    try {
        tBottom = state.at(tBottomName);
    } catch (const std::out_of_range& e) {
        Logged::info("No " + tBottomName
            + " field in restart file. Setting it to the melting point of ice.\n");
        tBottom = seaIceTf;
    }
}

/*
 * Use of double instead of FloatType is intentional here!
 * Some of these computations have caused problems in float, causing rapid ice growth in summer for
 * longer simulations.
 * The performance cost is only marginal. We still benefit from the effective
 * bandwith increase when the fields are stored as float and the thermodynamics are cheap in terms
 * of compute.
 */
KERNEL_IMPL_FUNCTION static void calculateTemps(double& tSurf, double& tUppr, double& tLowr,
    double& mSurf, const double cice, const double dQia_dt, const double hice, const double hsnow,
    const double penSw, const double qia, const double tf, double dt, const double cVol,
    const double seaIceTf, const double kappa_s)
{
    /*
     * Internal temperatures
     * Following Winton (2000) we start by calculating the temperature - first T1, then Tsurf, and
     * finally T2 Numers in parentheses refer to equations in the paper
     */

    const double hi = hice / cice;
    const double hs = hsnow / cice;
    const double tBase = tf; // Freezing point of seawater with the local salinity
    const double tMelt = (hs > 0) ? 0 : seaIceTf; // Melting point at the surface

    // First some coefficients based on temperatures from the previous time step
    double k12 = 4 * Ice::kappa * kappa_s / (kappa_s * hi + 4 * Ice::kappa * hs); // (5)
    double a = qia - tSurf * dQia_dt; // (7)
    double b = dQia_dt; // (8)
    double k32 = 2 * Ice::kappa / hi; // (10)

    double a1 = hi * cVol / (2 * dt) + k32 * (4 * dt * k32 + hi * cVol) / (6 * dt * k32 + hi * cVol)
        + b * k12 / (k12 + b); // (16)
    double b1 = -hi * (cVol * tUppr + Ice::Lf * Ice::rho * seaIceTf / tUppr) / (2 * dt) - penSw
        - k32 * (4 * dt * k32 * tBase + hi * cVol * tLowr) / (6 * dt * k32 + hi * cVol)
        + a * k12 / (k12 + b); // (17)
    double c1 = hi * Ice::Lf * Ice::rho * seaIceTf / (2 * dt); // (18)

    // Updated surface and mid-ice temperatures
    tUppr = -(b1 + Utils::sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1); // (21)
    tSurf = (k12 * tUppr - a) / (k12 + b); // (6)
    // Is the surface melting?
    if (tSurf > tMelt) {
        // Recalculate T1 and Tsurf if so
        tSurf = tMelt;
        // apply the change to the *1 parameters, rather than recalculating in full
        a1 += k12 - k12 * b / (k12 + b); // (19) simplified using +=
        b1 -= k12 * tSurf + a * k12 / (k12 + b); // (20) simplified using -=
        tUppr = -(b1 + Utils::sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1); // (21)

        // Surface melt
        mSurf = k12 * (tUppr - tSurf) - (a + b * tSurf); // (22)
    }

    // update the lower T based on the new value of the upper T (eq 15)
    tLowr = (2 * dt * k32 * (tUppr + 2 * tf) + hi * cVol * tLowr) / (6 * dt * k32 + hi * cVol);
}

void ThermoWinton::update(const TimestepTime& tst)
{
    // Advect ice temperatures
    IIceThermodynamics::update(tst);
    auto execSpace = DefaultExecutionSpace();
    auto& tBottom = tBottomAccessor.getAutoRW(execSpace);
    auto& tInternal = tInternalAccessor.getAutoRW(execSpace);
    FieldAdvection::advectField(tInternal, tst, IIceThermodynamics::minT, seaIceTf);
    FieldAdvection::advectField(tBottom, tst, IIceThermodynamics::minT, seaIceTf);

    // explicit execution space enables asynchronous execution
    auto& hsnow = hsnowAccessor.getAutoRW(execSpace);
    auto& botMelt = botMeltAccessor.getAutoRW(execSpace);
    auto& hice = hiceAccessor.getAutoRW(execSpace);
    auto& topMelt = topMeltAccessor.getAutoRW(execSpace);
    auto& qswBase = qswBaseAccessor.getAutoRW(execSpace);
    auto& snowMelt = snowMeltAccessor.getAutoRW(execSpace);
    auto& tsurf = tsurfAccessor.getAutoRW(execSpace);
    auto& deltaHi = deltaHiAccessor.getAutoRW(execSpace);
    auto& qow = qowAccessor.getAutoRW(execSpace);
    auto& snowToIce = snowToIceAccessor.getAutoRW(execSpace);
    auto& cice = ciceAccessor.getAutoRW(execSpace);
    const auto& dQia_dt = dQia_dtAccessor.getAutoRO(execSpace);
    const auto& penSw = penSwAccessor.getAutoRO(execSpace);
    const auto& qia = qiaAccessor.getAutoRO(execSpace);
    const auto& qio = qioAccessor.getAutoRO(execSpace);
    const auto& snowfall = snowfallAccessor.getAutoRO(execSpace);
    const auto& subl = sublAccessor.getAutoRO(execSpace);
    const auto& tf = tfAccessor.getAutoRO(execSpace);

    // static members can not be captured directly
    const double cVol = ThermoWinton::cVol;
    const bool doFlooding = ThermoWinton::doFlooding;
    const double seaIceTf = ThermoWinton::seaIceTf;
    const double kappa_s = ThermoWinton::kappa_s;
    const double cMin = IceMinima::c();
    const double hMin = IceMinima::h();
    const double dt = tst.step.seconds();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        static const double bulkLHFusionSnow = Water::Lf * Ice::rhoSnow;
        static const double bulkLHFusionIce = Water::Lf * Ice::rho;

        // Take care of too thin (or no) ice
        if (cice[i] <= cMin || hice[i] <= hMin) {
            if (cice[i] > 0) {
                deltaHi[i] = hice[i] / cice[i];
                botMelt[i] = topMelt[i] = 0.5 * deltaHi[i];
                snowMelt[i] = hsnow[i] / cice[i];
            } else {
                deltaHi[i] = 0.;
                botMelt[i] = topMelt[i] = 0.;
                snowMelt[i] = 0.;
            }

            snowToIce[i] = 0.;
            qswBase[i] = 0.;

            // Add to open water flux, since cice will be set to zero
            qow[i] += (hice[i] * bulkLHFusionIce + hsnow[i] * bulkLHFusionSnow) / dt;

            cice[i] = 0;
            hice[i] = 0;
            hsnow[i] = 0;

            tsurf[i] = seaIceTf;
            tInternal[i] = seaIceTf;
            tBottom[i] = seaIceTf;

            return;
        }

        double tSurf = tsurf[i]; // surface temperature
        double tUppr = tInternal[i]; // upper layer temperature
        double tLowr = tBottom[i]; // lower layer temperature
        double tBott = tf[i]; // freezing point of (local) seawater

        double hi = hice[i] / cice[i];
        double hs = hsnow[i] / cice[i];
        const double oldHi = hi; // Ice thickness at the start of the timestep

        double surfMelt = 0; // surface melting mass loss
        // Calculate temperatures by solving the heat conduction equation
        calculateTemps(tSurf, tUppr, tLowr, surfMelt, cice[i], dQia_dt[i], hice[i], hsnow[i],
            penSw[i], qia[i], tf[i], dt, cVol, seaIceTf, kappa_s);

        // The ratio of ΔΗ_f T_f / c_p,ice is used a lot. Units are K²
        const static double dHfTf_cp = Water::Lf * seaIceTf / Ice::cp;

        // Thickness changes
        // ice
        double h1 = hi / 2;
        double h2 = hi / 2;
        // Eqs. (1) and (25) - but I 've multiplied them with \rho_i (hence cVol), because it's
        // missing in the paper
        double e1 = cVol * (tUppr - seaIceTf) - bulkLHFusionIce * (1 - seaIceTf / tUppr);
        double e2 = cVol * (tLowr - seaIceTf) - bulkLHFusionIce;

        // snow
        hs += snowfall[i] / Ice::rhoSnow * dt;
        //    double accumulatedSnowThickness = snowfall[i] / Ice::rhoSnow * dt;

        // sublimation
        // 4 cases
        const double& subli = subl[i];
        double deltaSnow = subli * dt / Ice::rhoSnow;
        double deltaIce1 = (deltaSnow - hs) * Ice::rhoSnow / Ice::rho;
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
            // double oceanEvapError = (deltaIce2 - h2) * Ice::rho / Water::rhoOcean;
            // TODO: log the error
            h2 = 0;
            h1 = 0;
            hs = 0;
        }
        // Sublimated ice counts as top melt
        topMelt[i] = Utils::max(0.0, h1 + h2 - hi); // (23)

        // Bottom melt/freezing
        double meltBottom = (qio[i] - 4 * Ice::kappa * (tBott - tLowr) / hi) * dt;
        snowMelt[i] = 0;
        if (meltBottom <= 0.0) {
            // Freezing
            // Eq. (25) - but I 've multiplied them with \rho_i (hence cVol), because it's
            // missing in the paper
            double eBot = cVol * (tBott - seaIceTf) - bulkLHFusionIce;
            deltaIce2 = meltBottom / eBot; // (24)
            tLowr = (deltaIce2 * tBott + h2 * tLowr) / (deltaIce2 + h2); // (26)
            h2 += deltaIce2;
        } else {
            // Melting
            // Eqs. (31)-(32) with added division with \rho_i (and \rho_s for 32)
            deltaIce2 = -Utils::min(-meltBottom / e2, h2);
            deltaIce1 = -Utils::min(Utils::max(-(meltBottom + e2 * h2) / e1, 0.0), h1);
            snowMelt[i] = -Utils::min(
                Utils::max((meltBottom + e2 * h2 + e1 * h1) / bulkLHFusionSnow, 0.0), hs);

            // If everything melts we need to put heat back into the ocean
            if (h2 + h1 + hs - deltaIce2 - deltaIce1 - snowMelt[i] <= 0.0) {
                // (34) - with added multiplication of rhoi and rhos and division with dt
                qow[i] -= cice[i]
                    * Utils::max(meltBottom - bulkLHFusionSnow * hs + e1 * h1 + e2 * h2, 0.0) / dt;
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
        snowMelt[i] -= Utils::min(surfMelt * dt / bulkLHFusionSnow, hs);
        deltaIce1 = -Utils::min(Utils::max(-(surfMelt * dt - bulkLHFusionSnow * hs) / e1, 0.0), h1);
        deltaIce2 = -Utils::min(
            Utils::max(-(surfMelt * dt - bulkLHFusionSnow * hs + e1 * h1) / e2, 0.0), h2);

        // If everything melts we need to put heat back into the ocean
        // Eq (30) - with multiplication of rhoi and rhos and division with dt
        if (h2 + h1 + hs - deltaIce2 - deltaIce1 - snowMelt[i] <= 0.0) {
            qow[i] -= cice[i]
                * Utils::max(surfMelt * dt - bulkLHFusionSnow * hs + e1 * h1 + e2 * h2, 0.0) / dt;
        }

        hs += snowMelt[i];
        h1 += deltaIce1;
        h2 += deltaIce2;
        topMelt[i] += deltaIce1 + deltaIce2;

        // Snow to ice conversion
        double freeboard
            = (hi * (Water::rhoOcean - Ice::rho) - hs * Ice::rhoSnow) / Water::rhoOcean;
        if (doFlooding && freeboard < 0.0) {
            hs += Utils::min(freeboard * Ice::rho / Ice::rhoSnow, 0.0); // (35) using +=
            deltaIce1 = Utils::max(-freeboard, 0.0); // (36)
            double f1 = 1 - deltaIce1 / (deltaIce1 + h1); // Fraction of new ice in the upper layer
            double tBar = f1 * (tUppr + dHfTf_cp / tUppr) + (1 - f1) * seaIceTf; // (39)
            tUppr = (tBar - Utils::sqrt(tBar * tBar - 4 * dHfTf_cp)) / 2; // (38)
            h1 += deltaIce1;
            snowToIce[i] += deltaIce1;
        }

        // Add up the half-layer thicknesses
        hi = h1 + h2;
        // Adjust the temperatures to evenly divide the ice
        if (h2 > h1) {
            // Lower layer ice is added to the upper layer
            double f1 = h1 / hi * 2;
            double tBar = f1 * (tUppr + dHfTf_cp / tUppr) + (1 - f1) * tLowr; // (39)
            // The upper layer temperature changes
            tUppr = (tBar - Utils::sqrt(tBar * tBar - 4 * dHfTf_cp)) / 2; // (38)
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
        deltaHi[i] = hi - oldHi;

        // Remove very small ice thickness
        // The fluxes are already sorted
        if (hi < hMin) {
            if (deltaHi[i] < 0) {
                topMelt[i] *= oldHi / deltaHi[i];
                botMelt[i] *= oldHi / deltaHi[i];
            }
            snowToIce[i] = 0;

            deltaHi[i] = -oldHi;
            cice[i] = 0;

            tSurf = seaIceTf;
            tUppr = seaIceTf;
            tLowr = seaIceTf;
        }

        tsurf[i] = tSurf;
        tInternal[i] = tUppr;
        tBottom[i] = tLowr;

        hice[i] = hi * cice[i];
        hsnow[i] = hs * cice[i];
    });
}

} /* namespace Nextsim */
