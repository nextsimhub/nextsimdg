/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
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
    : hiceDG(getStore())
    , ciceDG(getStore())
    , hsnowDG(getStore())
    , qow(getStore())
    , deltaHi(getStore())
{
}

void IceGrowth::setData(const ModelState::DataMap& ms)
{
    iVertical->setData(ms);
    iLateral->setData(ms);
    iHealing->setData(ms);
}

ModelState IceGrowth::getStateDiagnostic() const
{
    ModelState state = iLateral->getStateDiagnostic();
    state.merge(iVertical->getStateDiagnostic());
    state.merge(iHealing->getStateDiagnostic());

    return state;
}

ModelState IceGrowth::getStatePrognostic() const
{
    ModelState state;
    // Merge in other states here
    state.merge(iLateral->getStatePrognostic());
    state.merge(iVertical->getStatePrognostic());
    state.merge(iHealing->getStatePrognostic());

    return state;
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
    // The snowMelt array is not currently filled with data, but it used elsewhere
    // FIXME calculate a true value for snowMelt

    if (doThermo) {
        iVertical->update(tsTime);
        iLateral->update(tsTime);
    }

    // Damage always heals, even if there's no active thermo
    // TODO: This should only be called for brittle rheologies
    iHealing->update(tsTime);
}
} /* namespace Nextsim */
