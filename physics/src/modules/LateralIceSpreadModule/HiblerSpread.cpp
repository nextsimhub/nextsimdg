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
KERNEL_IMPL_FUNCTION FloatType freeze(const FloatType newIce, const FloatType h0)
{
    return newIce / h0;
}

/*!
 * Updates the lateral melting of ice for the timestep.
 *
 * @param deltaHi The change in ice thickness this timestep.
 * @param cice The ice concentration.
 * @param hice The ice-average ice thickness.
 */
KERNEL_IMPL_FUNCTION FloatType melt(
    const FloatType deltaHi, const FloatType cice, const FloatType hice)
{
    /* We only decrease the concentration if the ice is melting, if the ice cover is not 100%, and
     * if there's ice there in the first place. */
    if (deltaHi < 0._ft && 0._ft < cice && cice < 1._ft) {
        const FloatType hiOld = hice / cice - deltaHi;
        return deltaHi * cice / (2._ft * hiOld);
    } else {
        return 0._ft;
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
    const FloatType h0Local = h0;
    const FloatType dt = tstep.step.seconds();
    const FloatType cMin = IceMinima::c();
    const FloatType hMin = IceMinima::h();

    overElementsAuto(OVER_ELEMENTS_LAMBDA(const ElementIndex i) {
        // newIceFormation
        // Flux cooling the ocean from open water
        // TODO Add assimilation fluxes here
        const FloatType coolingFlux = qow[i];
        // Temperature change of the mixed layer during this timestep
        const FloatType deltaTml = -coolingFlux / mixedLayerBulkHeatCapacity[i] * dt;
        // Initial temperature
        const FloatType t0 = sst[i];
        // Freezing point temperature
        const FloatType tf0 = tf[i];

        // deal with cooling below the freezing point
        if (const FloatType t1 = t0 + deltaTml; t1 < tf0) {
            // Heat lost cooling the mixed layer to freezing point
            const FloatType sensibleFlux = (tf0 - t0) / deltaTml * coolingFlux;
            // Any heat beyond that is latent heat forming new ice
            const FloatType latentFlux = coolingFlux - sensibleFlux;

            qow[i] = sensibleFlux;
            newIce[i] = latentFlux * dt * (1 - cIce[i]) / (Ice::Lf * Ice::rho);
        } else {
            newIce[i] = 0;
        }

        // lateralIceSpread
        const FloatType deltaCMelt = melt(deltaHi[i], cIce[i], hIce[i]);
        const FloatType deltaCFreeze = freeze(newIce[i], h0Local);
        deltaCIce[i] = deltaCFreeze + deltaCMelt;

        /* Snow is lost if the concentration decreases, and energy is returned to the ocean. We
         * reduce the snow volume by a "slice" of snow with the dimensions hs * deltaCIce.
         */
        const FloatType hs = hSnow[i] / std::max(cMin, cIce[i]);
        qow[i] -= deltaCMelt * hs * Water::Lf * Ice::rhoSnow / dt;
        hSnow[i] += hs * deltaCMelt;

        // Update ice concentration and volume
        cIce[i] += deltaCIce[i];
        hIce[i] += newIce[i];

        // applyLimits
        if (cIce[i] < cMin || hIce[i] < hMin) {
            qow[i] += Water::Lf * (hIce[i] * Ice::rho + hSnow[i] * Ice::rhoSnow) / dt;
            deltaCIce[i] = -cIce[i];
            deltaHice[i] = -hIce[i];
            deltaSmelt[i] = -hSnow[i];

            hIce[i] = 0;
            cIce[i] = 0;
            hSnow[i] = 0;
            newIce[i] = 0;
        }
    });
    timer.stop();
}
}
