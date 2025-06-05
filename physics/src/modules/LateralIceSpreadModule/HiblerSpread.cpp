/*!
 * @file HiblerSpread.cpp
 *
 * @date 04 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/HiblerSpread.hpp"

#include "include/IceMinima.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double HiblerSpread::h0 = 0;
double HiblerSpread::phiM = 0;

static const double h0Default = 0.25;
static const double phimDefault = 0.5;

static const std::map<int, std::string> keyMap = {
    { HiblerSpread::H0_KEY, "Hibler.h0" },
    { HiblerSpread::PHIM_KEY, "Hibler.phiM" },
};

void HiblerSpread::configure()
{
    h0 = Configured::getConfiguration(keyMap.at(H0_KEY), h0Default);
    phiM = Configured::getConfiguration(keyMap.at(PHIM_KEY), phimDefault);
}

ConfigMap HiblerSpread::getConfiguration() const
{
    return {
        { keyMap.at(H0_KEY), h0 },
        { keyMap.at(PHIM_KEY), phiM },
    };
}

ModelState HiblerSpread::getStateDiagnostic() const { return { {}, getConfiguration() }; }

HiblerSpread::HelpMap& HiblerSpread::getHelpText(HelpMap& map, bool getAll)
{
    map["HiblerSpread"] = {
        { keyMap.at(H0_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(h0Default), "m", "The thickness of newly frozen ice." },
        { keyMap.at(PHIM_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(phimDefault), "", "Power-law exponent for melting ice." },
    };
    return map;
}
HiblerSpread::HelpMap& HiblerSpread::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

void HiblerSpread::freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
    double newIce, double& cice, double& qow, double& deltaCfreeze)
{
    static const double ooh0 = 1. / h0;
    deltaCfreeze = newIce * ooh0;
}

void HiblerSpread::melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
    double& cice, double& qow, double& deltaCmelt)
{
    if (cice < 1) {
        deltaCmelt = deltaHi * cice * phiM / hice;
    }
}

void HiblerSpread::newIceFormation(size_t i, const TimestepTime& tst)
{
    // Flux cooling the ocean from open water
    // TODO Add assimilation fluxes here
    double coolingFlux = qow[i];
    // Temperature change of the mixed layer during this timestep
    double deltaTml = -coolingFlux / mixedLayerBulkHeatCapacity[i] * tst.step;
    // Initial temperature
    double t0 = sst[i];
    // Freezing point temperature
    double tf0 = tf[i];
    // Final temperature
    double t1 = t0 + deltaTml;

    // deal with cooling below the freezing point
    if (t1 < tf0) {
        // Heat lost cooling the mixed layer to freezing point
        double sensibleFlux = (tf0 - t0) / deltaTml * coolingFlux;
        // Any heat beyond that is latent heat forming new ice
        double latentFlux = coolingFlux - sensibleFlux;

        qow[i] = sensibleFlux;
        newice[i] = latentFlux * tst.step * (1 - cice[i]) / (Ice::Lf * Ice::rho);
    } else {
        newice[i] = 0;
    }
}

// Update thickness with concentration
double HiblerSpread::updateThickness(double& thick, double newConc, double deltaC, double deltaV)
{
    return thick += (deltaV - thick * deltaC) / newConc;
}

void HiblerSpread::lateralIceSpread(size_t i, const TimestepTime& tstep)
{
    double deltaCMelt = 0;
    double deltaCFreeze = 0;
    freeze(tstep, hice[i], hsnow[i], deltaHi[i], newice[i], cice[i], qow[i], deltaCFreeze);
    if (deltaHi[i] < 0) {
        // Note that the cell-averaged hice0 is converted to an ice averaged value
        melt(tstep, hice0[i], hsnow[i], deltaHi[i], cice[i], qow[i], deltaCMelt);
    }
    deltaCIce[i] = deltaCFreeze + deltaCMelt;
    cice[i] = (hice[i] > 0 || newice[i] > 0) ? cice[i] + deltaCIce[i] : 0;
    if (cice[i] >= IceMinima::c()) {
        // The updated ice thickness must conserve volume
        updateThickness(hice[i], cice[i], deltaCIce[i], newice[i]);
        if (deltaCIce[i] < 0) {
            // Snow is lost if the concentration decreases, and energy is returned to the ocean
            qow[i] -= deltaCIce[i] * hsnow[i] * Water::Lf * Ice::rhoSnow / tstep.step;
        } else {
            // Update snow thickness. Currently, no new snow is implemented
            updateThickness(hsnow[i], cice[i], deltaCIce[i], 0);
        }
    }
}

void HiblerSpread::applyLimits(size_t i, const TimestepTime& tstep)
{
    if (cice[i] < IceMinima::c() || hice[i] < IceMinima::h()) {
        qow[i] += cice[i] * Water::Lf * (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / tstep.step;
        hice[i] = 0;
        cice[i] = 0;
        hsnow[i] = 0;
    }
}
}
