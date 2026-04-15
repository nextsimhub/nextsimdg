/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IDynamics.hpp"
#include "include/NextsimModule.hpp"

#include "PDTestDynamics.hpp"

#include <string>

namespace Module {
const std::string PDTESTDYNAMICS = "Nextsim::PDTestDynamics";

template <> const Module<Nextsim::IDynamics>::Map& Module<Nextsim::IDynamics>::functionMap()
{
    static const Map theMap = {
        { PDTESTDYNAMICS, newImpl<Nextsim::IDynamics, Nextsim::PDTestDynamics> },
    };
    return theMap;
}

// Needed so that the finalizer works correctly because getGenerationFunction() looks up the
// generator function under this name. For some reason the PDTestDynamics module is only registered
// in a build without debug symbols, so this issue only appears in a release build.
template <> std::string Module<Nextsim::IDynamics>::getDefaultImplementationName()
{
    return PDTESTDYNAMICS;
}

template <> Module<Nextsim::IDynamics>::Fn& Module<Nextsim::IDynamics>::getGenerationFunction()
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
