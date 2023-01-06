/*!
 * @file DynamicsModule.cpp
 *
 * @date 6 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DynamicsModule.hpp"

#include "include/DummyDynamics.hpp"

#include <string>

namespace Module {
const std::string DUMMYDYNAMICS = "Nextsim::DummyDynamics";

template <>
Module<Nextsim::IDynamics>::map Module<Nextsim::IDynamics>::functionMap = {
    { DUMMYDYNAMICS, newImpl<Nextsim::IDynamics, Nextsim::DummyDynamics> },
};

template <>
Module<Nextsim::IDynamics>::fn Module<Nextsim::IDynamics>::spf = functionMap.at(DUMMYDYNAMICS);
template <>
std::unique_ptr<Nextsim::IDynamics> Module<Nextsim::IDynamics>::staticInstance
= std::move(newImpl<Nextsim::IDynamics, Nextsim::DummyDynamics>());

template <>
std::string Module<Nextsim::IDynamics>::moduleName(){ return "Nextsim::IDynamics"; }

template<> HelpMap& getHelpRecursive<Nextsim::IDynamics>(HelpMap& map, bool getAll)
{
    const std::string& pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IDynamics>::moduleName(), ConfigType::MODULE,
        { DUMMYDYNAMICS }, DUMMYDYNAMICS, "",
        "MODULE DESCRIPTION HERE" });
    return map;
}
template <>
Nextsim::IDynamics& getImplementation<Nextsim::IDynamics>()
{
    return getImplTemplate<Nextsim::IDynamics, DynamicsModule>();
}
template <>
void setImplementation<Nextsim::IDynamics>(const std::string& implName)
{
    setImplTemplate<DynamicsModule>(implName);
}
template <>
std::unique_ptr<Nextsim::IDynamics> getInstance()
{
    return getInstTemplate<Nextsim::IDynamics, DynamicsModule>();
}
DynamicsModule::Constructor DynamicsModule::ctor;
DynamicsModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IDynamics, DynamicsModule>();
}

} /* namespace Module */
