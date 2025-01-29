/*!
 * @file IceGrowth.cpp
 *
 * @date 29 Jan 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/IceGrowth.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

static const std::map<int, std::string> keyMap = {
    { IceGrowth::ICE_THERMODYNAMICS_KEY, "IceThermodynamicsModel" },
    { IceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel" },
    { IceGrowth::MINC_KEY, "nextsim_thermo.min_conc" },
    { IceGrowth::MINH_KEY, "nextsim_thermo.min_thick" },
    { IceGrowth::USE_THERMO_KEY, "nextsim_thermo.use_thermo_forcing" },
};

IceGrowth::IceGrowth()
    : deltaCFreeze(ModelArray::Type::H)
    , deltaCIce(ModelArray::Type::H)
    , deltaCMelt(ModelArray::Type::H)
    , newice(ModelArray::Type::H)
    , cice(getStore())
    , deltaHi(getStore())
    , hice(getStore())
    , hsnow(getStore())
    , mixedLayerBulkHeatCapacity(getStore())
    , sst(getStore())
    , tf(getStore())
    , qow(getStore())
{
    registerModule();
    getStore().registerArray(Shared::NEW_ICE, &newice, RW);
    getStore().registerArray(Shared::DELTA_CICE, &deltaCIce, RW);
}

void IceGrowth::setData(const ModelState::DataMap& ms)
{
    newice.resize();
    deltaCFreeze.resize();
    deltaCMelt.resize();
    deltaCIce.resize();

    iVertical->setData(ms);
    iLateral->setData(ms);
    iHealing->setData(ms);
}

ModelState IceGrowth::getState() const { return { {}, getConfiguration() }; }

ModelState IceGrowth::getStateRecursive(const OutputSpec& os) const
{
    ModelState state(getState());
    // Merge in other states here
    state.merge(iLateral->getStateRecursive(os));
    state.merge(iVertical->getStateRecursive(os));
    state.merge(iHealing->getStateRecursive(os));

    return os ? state : ModelState();
}

IceGrowth::HelpMap& IceGrowth::getHelpText(HelpMap& map, bool getAll)
{
    map["IceGrowth"] = {
        { keyMap.at(MINC_KEY), ConfigType::NUMERIC, { "0", "1" },
            ConfigurationHelp::toString(IceMinima::cMinDefault), "",
            "Minimum allowed ice concentration." },
        { keyMap.at(MINH_KEY), ConfigType::NUMERIC, { "0", "∞" },
            ConfigurationHelp::toString(IceMinima::hMinDefault), "m",
            "Minimum allowed ice thickness." },
        { keyMap.at(USE_THERMO_KEY), ConfigType::BOOLEAN, { "true", "false" }, "true", "",
            "Perform ice physics calculations as part of the timestep." },
    };
    return map;
}
IceGrowth::HelpMap& IceGrowth::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    Module::getHelpRecursive<IIceThermodynamics>(map, getAll);
    Module::getHelpRecursive<ILateralIceSpread>(map, getAll);
    Module::getHelpRecursive<IDamageHealing>(map, getAll);
    return map;
}

void IceGrowth::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceThermodynamics>);
    Finalizer::registerUnique(Module::finalize<ILateralIceSpread>);
    Finalizer::registerUnique(Module::finalize<IDamageHealing>);

    // Configure whether we actually do anything here
    doThermo = Configured::getConfiguration(keyMap.at(USE_THERMO_KEY), true);
    // Configure constants
    IceMinima::cMin = Configured::getConfiguration(keyMap.at(MINC_KEY), IceMinima::cMinDefault);
    IceMinima::hMin = Configured::getConfiguration(keyMap.at(MINH_KEY), IceMinima::hMinDefault);

    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IIceThermodynamics>());
    iLateral = std::move(Module::getInstance<ILateralIceSpread>());
    iHealing = std::move(Module::getInstance<IDamageHealing>());
    tryConfigure(*iVertical);
    tryConfigure(*iLateral);
    tryConfigure(*iHealing);
}

ConfigMap IceGrowth::getConfiguration() const
{
    return {
        { keyMap.at(MINC_KEY), IceMinima::cMin },
        { keyMap.at(MINH_KEY), IceMinima::hMin },
    };
}

void IceGrowth::update(const TimestepTime& tsTime)
{
    if (doThermo) {
        iVertical->update(tsTime);
        // new ice formation
        overElements(std::bind(&IceGrowth::updateWrapper, this, std::placeholders::_1,
                         std::placeholders::_2),
            tsTime);
    }

    // Damage always heals, even if there's no active thermo
    // TODO: This should only be called for brittle rheologies
    iHealing->update(tsTime);

    // Apply sensible limits
    overElements(
        std::bind(&IceGrowth::applyLimits, this, std::placeholders::_1, std::placeholders::_2),
        tsTime);
}

void IceGrowth::newIceFormation(size_t i, const TimestepTime& tst)
{
    // Flux cooling the ocean from open water
    // TODO Add assimilation fluxes here
    const double coolingFlux = qow[i];
    // Temperature change of the mixed layer during this timestep
    const double deltaTml = -coolingFlux / mixedLayerBulkHeatCapacity[i] * tst.step;
    // Initial temperature
    const double t0 = sst[i];
    // Freezing point temperature
    const double tf0 = tf[i];
    // Final temperature
    const double t1 = t0 + deltaTml;

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

void IceGrowth::lateralIceSpread(size_t i, const TimestepTime& tstep)
{
    deltaCMelt[i] = 0;
    deltaCFreeze[i] = 0;
    const double hsOld = (cice[i] > 0) ? hsnow[i] / cice[i] : 0;
    iLateral->freeze(
        tstep, hice[i], hsnow[i], deltaHi[i], newice[i], cice[i], qow[i], deltaCFreeze[i]);
    if (deltaHi[i] < 0) {
        iLateral->melt(tstep, hice[i], hsnow[i], deltaHi[i], cice[i], qow[i], deltaCMelt[i]);
    }
    deltaCIce[i] = deltaCFreeze[i] + deltaCMelt[i];
    cice[i] = (hice[i] > 0 || newice[i] > 0) ? cice[i] + deltaCIce[i] : 0;
    if (cice[i] >= IceMinima::c()) {
        // The updated ice volume
        hice[i] += newice[i];
        if (deltaCIce[i] < 0) {
            // Snow is lost if the concentration decreases, and energy is returned to the ocean
            qow[i] -= deltaCIce[i] * hsnow[i] * Water::Lf * Ice::rhoSnow / tstep.step;
            // Conserve the snow (slab) thickness
            hsnow[i] = hsOld * cice[i];
        }
    }
}

void IceGrowth::applyLimits(size_t i, const TimestepTime& tstep)
{
    if (cice[i] < IceMinima::c() || hice[i] < IceMinima::h()) {
        qow[i] += cice[i] * Water::Lf * (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / tstep.step;
        hice[i] = 0;
        cice[i] = 0;
        hsnow[i] = 0;
    }
}
} /* namespace Nextsim */
