/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/HiblerSpread.hpp"
#include "include/KernelAlternatives.hpp"
#include "kokkos/include/KokkosTimer.hpp"

#include "include/IceMinima.hpp"
#include "include/constants.hpp"

namespace Nextsim {

FloatType HiblerSpread::h0 = 0;

static constexpr FloatType h0Default = 0.25;

static const std::map<int, std::string> keyMap = {
    { HiblerSpread::H0_KEY, "Hibler.h0" },
};

void HiblerSpread::configure() { h0 = Configured::getConfiguration(keyMap.at(H0_KEY), h0Default); }

ConfigMap HiblerSpread::getConfiguration() const
{
    return {
        { keyMap.at(H0_KEY), h0 },
    };
}

ModelState HiblerSpread::getStateDiagnostic() const { return { {}, getConfiguration() }; }

HiblerSpread::HelpMap& HiblerSpread::getHelpText(HelpMap& map, bool getAll)
{
    map["HiblerSpread"] = {
        { keyMap.at(H0_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(h0Default), "m", "The thickness of newly frozen ice." },
    };
    return map;
}
HiblerSpread::HelpMap& HiblerSpread::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

/*!
 * Updates the freezing of open water for the timestep.
 *
 * @param newIce The positive change in ice thickness this timestep.
 * @param h0 The demarcation thickness between thick and thin ice
 */
KERNEL_IMPL_FUNCTION FloatType freeze(const double newIce, const double h0) { return newIce / h0; }

/*!
 * Updates the lateral melting of ice for the timestep.
 *
 * @param deltaHi The change in ice thickness this timestep.
 * @param cice The ice concentration.
 * @param hice The ice-average ice thickness.
 */
KERNEL_IMPL_FUNCTION FloatType melt(const double deltaHi, const double cice, const double hice)
{
    /* We only decrease the concentration if the ice is melting, if the ice cover is not 100%, and
     * if there's ice there in the first place. */
    if (deltaHi < 0. && 0. < cice && cice < 1.) {
        const double hiOld = (hice - deltaHi) / cice;
        return deltaHi / (2. * hiOld);
    } else {
        return 0.;
    }
}

void HiblerSpread::update(const TimestepTime& tstep)
{
    static KokkosTimer<true> timer("HiblerSpread");
    timer.start();

    const auto execSpace = DefaultExecutionSpace();
    auto& hSnow = hsnowAccessor.getAutoRW(execSpace);
    auto& qow = qowAccessor.getAutoRW(execSpace);
    auto& cIce = ciceAccessor.getAutoRW(execSpace);
    auto& newIce = newiceAccessor.getAutoRW(execSpace);
    auto& hIce = hiceAccessor.getAutoRW(execSpace);
    auto& deltaCIce = deltaCIceAccessor.getAutoRW(execSpace);
    auto& deltaHice = deltaHiceAccessor.getAutoRW();
    auto& deltaSmelt = deltaSmeltAccessor.getAutoRW();

    const auto& mixedLayerBulkHeatCapacity
        = mixedLayerBulkHeatCapacityAccessor.getAutoRO(execSpace);
    const auto& deltaHi = deltaHiAccessor.getAutoRO(execSpace);
    const auto& tf = tfAccessor.getAutoRO(execSpace);
    const auto& sst = sstAccessor.getAutoRO(execSpace);

    // static members can not be captured directly
    const double h0Local = h0;
    const double dt = tstep.step.seconds();
    const double cMin = IceMinima::c();
    const double hMin = IceMinima::h();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // newIceFormation
        // Flux cooling the ocean from open water
        // TODO Add assimilation fluxes here
        double coolingFlux = qow[i];
        const double mixedLayerBulkHeatCapacityD = mixedLayerBulkHeatCapacity[i];
        // Temperature change of the mixed layer during this timestep
        const double deltaTml = -coolingFlux / mixedLayerBulkHeatCapacityD * dt;
        // Initial temperature
        const double t0 = sst[i];
        // Freezing point temperature
        const double tf0 = tf[i];

        double cice = cIce[i];
        double newIceD;

        // deal with cooling below the freezing point
        if (const double t1 = t0 + deltaTml; t1 < tf0) {
            // Heat lost cooling the mixed layer to freezing point.
            // Guard the 0/0 when there is no cooling flux.
            const double sensibleFlux = deltaTml ? (tf0 - t0) / deltaTml * coolingFlux : 0.0;
            // Any heat beyond that is latent heat forming new ice
            const double latentFlux = coolingFlux - sensibleFlux;

            coolingFlux = sensibleFlux;
            newIceD = latentFlux * dt * (1 - cice) / (Ice::Lf * Ice::rho);
        } else {
            newIceD = 0.0;
        }

        // lateralIceSpread
        double hice = hIce[i];
        const double deltaCMelt = melt(deltaHi[i], cice, hice);
        const double deltaCFreeze = freeze(newIceD, h0Local);
        const double deltaCIceD = deltaCFreeze + deltaCMelt;

        /* Snow is lost if the concentration decreases, and energy is returned to the ocean. We
         * reduce the snow volume by a "slice" of snow with the dimensions hs * deltaCIce.
         */
        double hSnowD = hSnow[i];
        const double hs = hSnowD / Utils::max(cMin, cice);
        coolingFlux -= deltaCMelt * hs * Water::Lf * Ice::rhoSnow / dt;
        hSnowD += hs * deltaCMelt;

        // Update ice concentration and volume
        cice += deltaCIceD;
        hice += newIceD;

        // applyLimits
        if (cice < cMin || hice < hMin) {
            coolingFlux += Water::Lf * (hice * Ice::rho + hSnowD * Ice::rhoSnow) / dt;
            deltaCIce[i] = -cice;
            deltaHice[i] = -hice;
            deltaSmelt[i] = -hSnowD;

            hIce[i] = 0;
            cIce[i] = 0;
            hSnow[i] = 0;
            newIce[i] = 0;
        } else {
            deltaCIce[i] = deltaCIceD;

            hIce[i] = hice;
            cIce[i] = cice;
            hSnow[i] = hSnowD;
            newIce[i] = newIceD;
        }
        qow[i] = coolingFlux;
    });
    timer.stop();
}
}
