/*!
 * @file DynamicsModuleForPDtest.cpp
 *
 * @date 5 May 2023
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

template <> std::string Module<Nextsim::IDynamics>::moduleName() { return "Nextsim::IDynamics"; }

template <> HelpMap& getHelpRecursive<Nextsim::IDynamics>(HelpMap& map, bool getAll) { return map; }
template <> Nextsim::IDynamics& getImplementation<Nextsim::IDynamics>()
{
    return getImplTemplate<Nextsim::IDynamics, DynamicsModule>();
}
template <> void setImplementation<Nextsim::IDynamics>(const std::string& implName)
{
    setImplTemplate<DynamicsModule>(implName);
}
template <> std::unique_ptr<Nextsim::IDynamics> getInstance()
{
    return getInstTemplate<Nextsim::IDynamics, DynamicsModule>();
}
DynamicsModule::Constructor DynamicsModule::ctor;
DynamicsModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IDynamics, DynamicsModule>();
}

} /* namespace Module */
