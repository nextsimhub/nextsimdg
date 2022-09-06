/*!
 * @file AtmosphereStateModule.cpp
 *
 * @date May 11, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/AtmosphereStateModule.hpp"

#include "include/DummyAtmosphereState.hpp"

#include <string>

namespace Module {
const std::string DUMMYATMOSPHERESTATE = "Nextsim::DummyAtmosphereState";

template <>
Module<Nextsim::AtmosphereState>::map Module<Nextsim::AtmosphereState>::functionMap = {
    { DUMMYATMOSPHERESTATE, newImpl<Nextsim::AtmosphereState, Nextsim::DummyAtmosphereState> },
};

template <>
Module<Nextsim::AtmosphereState>::fn Module<Nextsim::AtmosphereState>::spf
    = functionMap.at(DUMMYATMOSPHERESTATE);
template <>
std::unique_ptr<Nextsim::AtmosphereState> Module<Nextsim::AtmosphereState>::staticInstance
    = std::move(newImpl<Nextsim::AtmosphereState, Nextsim::DummyAtmosphereState>());

template <> std::string Module<Nextsim::AtmosphereState>::moduleName()
{
    return "Nextsim::AtmosphereState";
}

template <> HelpMap& getHelpRecursive<Nextsim::AtmosphereState>(HelpMap& map, bool getAll)
{
    const std::string pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({ pfx + "." + Module<Nextsim::AtmosphereState>::moduleName(),
        ConfigType::MODULE, { DUMMYATMOSPHERESTATE }, DUMMYATMOSPHERESTATE, "",
        "The module selecting how the state of the atmosphere is obtained." });
    return map;
}
template <> Nextsim::AtmosphereState& getImplementation<Nextsim::AtmosphereState>()
{
    return getImplTemplate<Nextsim::AtmosphereState, AtmosphereStateModule>();
}
template <> void setImplementation<Nextsim::AtmosphereState>(const std::string& implName)
{
    setImplTemplate<AtmosphereStateModule>(implName);
}
template <> std::unique_ptr<Nextsim::AtmosphereState> getInstance()
{
    return getInstTemplate<Nextsim::AtmosphereState, AtmosphereStateModule>();
}
AtmosphereStateModule::Constructor AtmosphereStateModule::ctor;
AtmosphereStateModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::AtmosphereState, AtmosphereStateModule>();
}

} /* namespace Module */
