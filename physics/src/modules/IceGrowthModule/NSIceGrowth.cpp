/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#include "include/NSIceGrowth.hpp"

#include "include/Finalizer.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

static const std::map<int, std::string> keyMap = {
    { NSIceGrowth::ICE_THERMODYNAMICS_KEY, "IceThermodynamicsModel" },
    { NSIceGrowth::LATERAL_GROWTH_KEY, "LateralIceModel" },
    { NSIceGrowth::MINC_KEY, "nextsim_thermo.min_conc" },
    { NSIceGrowth::MINH_KEY, "nextsim_thermo.min_thick" },
    { NSIceGrowth::USE_THERMO_KEY, "nextsim_thermo.use_thermo_forcing" },
};

NSIceGrowth::NSIceGrowth()
    : hice(getStore())
    , cice(getStore())
    , hsnow(getStore())
    , qow(getStore())
    , deltaHi(getStore())
{
}

void NSIceGrowth::setData(const ModelState::DataMap& ms)
{
    iVertical->setData(ms);
    iLateral->setData(ms);
    iHealing->setData(ms);
}

ModelState NSIceGrowth::getStateDiagnostic() const
{
    ModelState state = iLateral->getStateDiagnostic();
    state.merge(iVertical->getStateDiagnostic());
    state.merge(iHealing->getStateDiagnostic());

    return state;
}

ModelState NSIceGrowth::getStatePrognostic() const
{
    ModelState state;
    // Merge in other states here
    state.merge(iLateral->getStatePrognostic());
    state.merge(iVertical->getStatePrognostic());
    state.merge(iHealing->getStatePrognostic());

    return state;
}

NSIceGrowth::HelpMap& NSIceGrowth::getHelpText(HelpMap& map, bool getAll)
{
    map["NSIceGrowth"] = {
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
NSIceGrowth::HelpMap& NSIceGrowth::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    Module::getHelpRecursive<IIceThermodynamics>(map, getAll);
    Module::getHelpRecursive<ILateralIceSpread>(map, getAll);
    Module::getHelpRecursive<IDamageHealing>(map, getAll);
    return map;
}

void NSIceGrowth::configure()
{
    Finalizer::registerUnique(Module::finalize<IIceThermodynamics>);
    Finalizer::registerUnique(Module::finalize<ILateralIceSpread>);
    Finalizer::registerUnique(Module::finalize<IDamageHealing>);

    // Configure whether we actually do anything here
    doThermo = Configured::getConfiguration(keyMap.at(USE_THERMO_KEY), true);
    // Configure constants. Negative values trigger the default values in IIceGrowth
    IIceGrowth::setCMin(Configured::getConfiguration(keyMap.at(MINC_KEY), IceMinima::cMinDefault));
    IIceGrowth::setHMin(Configured::getConfiguration(keyMap.at(MINH_KEY), IceMinima::hMinDefault));

    // Configure the vertical and lateral growth modules
    iVertical = std::move(Module::getInstance<IIceThermodynamics>());
    iLateral = std::move(Module::getInstance<ILateralIceSpread>());
    iHealing = std::move(Module::getInstance<IDamageHealing>());
    tryConfigure(*iVertical);
    tryConfigure(*iLateral);
    tryConfigure(*iHealing);
}

ConfigMap NSIceGrowth::getConfiguration() const
{
    return {
        { keyMap.at(MINC_KEY), IceMinima::c() },
        { keyMap.at(MINH_KEY), IceMinima::h() },
    };
}

void NSIceGrowth::update(const TimestepTime& tsTime)
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
