/*!
 * @file DynamicsModuleForPDtest.cpp
 *
 * @date 5 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IDynamics.hpp"
#include "include/NextsimModule.hpp"

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

template <> std::string Module<Nextsim::IDynamics>::moduleName() { return "Nextsim::IDynamics"; }

template <> HelpMap& getHelpRecursive<Nextsim::IDynamics>(HelpMap& map, bool getAll) { return map; }

} /* namespace Module */
