/*!
 * @file IceGrowth.cpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceGrowth.hpp"

#include "include/VerticalIceGrowthModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

double IceGrowth::minc;
double IceGrowth::minh;

template <>
const std::map<int, std::string> Configured<IceGrowth>::keyMap = {
    { IceGrowth::VERTICAL_GROWTH_KEY, "VerticalIceModel" },
    { IceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel" },
    { IceGrowth::MINC_KEY, "nextsim_thermo.min_conc" },
    { IceGrowth::MINH_KEY, "nextsim_thermo.min_thick" },

};

void IceGrowth::configure()
{
    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IVerticalIceGrowth>());
    iLateral = std::move(Module::getInstance<ILateralIceSpread>());

    // Configure constants
    minc = Configured::getConfiguration(keyMap.at(MINC_KEY), 1e-12);
    minh = Configured::getConfiguration(keyMap.at(MINH_KEY), 0.01);
}

void IceGrowth::update(const TimestepTime& tsTime)
{
    hice = hice0;
    cice = cice0;
    hsnow = hsnow0;

    iVertical->update(tsTime);

    // new ice formation
    overElements(
        std::bind(&IceGrowth::updateWrapper, this, std::placeholders::_1, std::placeholders::_2),
        tsTime);
}

void IceGrowth::newIceFormation(size_t i, const TimestepTime& tst)
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
        newice[0] = 0;
    }
}

// Update thickness with concentration
static double updateThickness(double& thick, double newConc, double deltaC, double deltaV)
{
    return thick += (deltaV - thick * deltaC) / newConc;
}

void IceGrowth::lateralIceSpread(size_t i, const TimestepTime& tstep)
{
    iLateral->freeze(
        tstep, hice[i], hsnow[i], deltaHi[i], newice[i], cice[i], qow[i], deltaCFreeze[i]);
    if (deltaHi[i] < 0) {

        iLateral->melt(tstep, hice[i], hsnow[i], deltaHi[i], cice[i], qow[i], deltaCMelt[i]);
    }
    double deltaC = deltaCFreeze[i] + deltaCMelt[i];
    cice[i] += deltaC;

    if (cice[i] >= minc) {
        // The updated ice thickness must conserve volume
        updateThickness(hice[i], cice[i], deltaC, newice[i]);
        if (deltaC < 0) {
            // Snow is lost if the concentration decreases, and energy is returned to the ocean
            qow[i] -= deltaC * hsnow[i] * Water::Lf * Ice::rhoSnow / tstep.step;
        } else {
            // Update snow thickness. Currently no new snow is implemented
            updateThickness(hsnow[i], cice[i], deltaC, 0);
        }
    }
}

void IceGrowth::applyLimits(size_t i, const TimestepTime& tstep)
{
    if (cice[i] < minc || hice[i] < minh) {
        qow[i] += cice[i] * Water::Lf * (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / tstep.step;
        hice[i] = 0;
        cice[i] = 0;
        hsnow[i] = 0;
    }
}
} /* namespace Nextsim */
