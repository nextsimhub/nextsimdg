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
const Module<Nextsim::IDynamics>::map& Module<Nextsim::IDynamics>::functionMap()
{
    static const map theMap = {
            { DUMMYDYNAMICS, newImpl<Nextsim::IDynamics, Nextsim::DummyDynamics> },
    };
    return theMap;
}

template <>
Module<Nextsim::IDynamics>::fn& Module<Nextsim::IDynamics>::getGenerationFunction()
{
    static fn thePtr = functionMap().at(DUMMYDYNAMICS);
    return thePtr;
}
//template <>
//std::unique_ptr<Nextsim::IDynamics> Module<Nextsim::IDynamics>::staticInstance
//    = std::move(newImpl<Nextsim::IDynamics, Nextsim::DummyDynamics>());

template <> std::string Module<Nextsim::IDynamics>::moduleName() { return "Nextsim::IDynamics"; }

template <> HelpMap& getHelpRecursive<Nextsim::IDynamics>(HelpMap& map, bool getAll) { return map; }
template <> Nextsim::IDynamics& getImplementation<Nextsim::IDynamics>()
{
    return Module<Nextsim::IDynamics>::getImplementation();
}
template <> void setImplementation<Nextsim::IDynamics>(const std::string& implName)
{
    Module<Nextsim::IDynamics>::setImplementation(implName);
}
template <> std::unique_ptr<Nextsim::IDynamics> getInstance()
{
    return Module<Nextsim::IDynamics>::getInstance();
}
//DynamicsModule::Constructor DynamicsModule::ctor;
//DynamicsModule::Constructor::Constructor()
//{
//    addToConfiguredModules<Nextsim::IDynamics, DynamicsModule>();
//}

} /* namespace Module */
