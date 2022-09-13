/*!
 * @file OceanStateModule.cpp
 *
 * @date May 11, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/OceanStateModule.hpp"

#include "include/ConfiguredOcean.hpp"
#include "include/DummyOceanState.hpp"

#include <string>

namespace Module {
const std::string DUMMYOCEANSTATE = "Nextsim::DummyOceanState";
const std::string CONFIGUREDOCEAN = "Nextsim::ConfiguredOcean";

template <>
Module<Nextsim::OceanState>::map Module<Nextsim::OceanState>::functionMap = {
    { DUMMYOCEANSTATE, newImpl<Nextsim::OceanState, Nextsim::DummyOceanState> },
    { CONFIGUREDOCEAN, newImpl<Nextsim::OceanState, Nextsim::ConfiguredOcean> },
};

template <>
Module<Nextsim::OceanState>::fn Module<Nextsim::OceanState>::spf = functionMap.at(DUMMYOCEANSTATE);
template <>
std::unique_ptr<Nextsim::OceanState> Module<Nextsim::OceanState>::staticInstance
    = std::move(newImpl<Nextsim::OceanState, Nextsim::DummyOceanState>());

template <> std::string Module<Nextsim::OceanState>::moduleName() { return "Nextsim::OceanState"; }

template <> HelpMap& getHelpRecursive<Nextsim::OceanState>(HelpMap& map, bool getAll)
{
    const std::string pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({ pfx + "." + Module<Nextsim::OceanState>::moduleName(), ConfigType::MODULE,
        { DUMMYOCEANSTATE }, DUMMYOCEANSTATE, "",
        "The module selecting how the state of the ocean is obtained." });
    return map;
}
template <> Nextsim::OceanState& getImplementation<Nextsim::OceanState>()
{
    return getImplTemplate<Nextsim::OceanState, OceanStateModule>();
}
template <> void setImplementation<Nextsim::OceanState>(const std::string& implName)
{
    setImplTemplate<OceanStateModule>(implName);
}
template <> std::unique_ptr<Nextsim::OceanState> getInstance()
{
    return getInstTemplate<Nextsim::OceanState, OceanStateModule>();
}
OceanStateModule::Constructor OceanStateModule::ctor;
OceanStateModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::OceanState, OceanStateModule>();
}

} /* namespace Module */
