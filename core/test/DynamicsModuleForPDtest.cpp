/*!
 * @file DynamicsModuleForPDtest.cpp
 *
 * @date 5 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IDynamics.hpp"
#include "include/NextsimModule.hpp"

#include "PDTestDynamics.hpp"

#include <string>

namespace Module {
const std::string PDTESTDYNAMICS = "Nextsim::PDTestDynamics";

template <>
const Module<Nextsim::IDynamics>::Map& Module<Nextsim::IDynamics>::functionMap()
{
    static const Map theMap = {
            { PDTESTDYNAMICS, newImpl<Nextsim::IDynamics, Nextsim::PDTestDynamics> },
    };
    return theMap;
}

template <>
Module<Nextsim::IDynamics>::Fn& Module<Nextsim::IDynamics>::getGenerationFunction()
{
    static Fn thePtr = functionMap().at(PDTESTDYNAMICS);
    return thePtr;
}

template <> std::string Module<Nextsim::IDynamics>::moduleName() { return "Nextsim::IDynamics"; }

template <> Nextsim::IDynamics& getImplementation<Nextsim::IDynamics>()
{
    return Module<Nextsim::IDynamics>::getImplementation();
}

template <> HelpMap& getHelpRecursive<Nextsim::IDynamics>(HelpMap& map, bool getAll) { return map; }

} /* namespace Module */
