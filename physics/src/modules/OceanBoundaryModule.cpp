/*!
 * @file OceanBoundaryModule.cpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/OceanBoundaryModule.hpp"

#include "include/ConstantOceanBoundary.hpp"

#include <string>

namespace Module {
const std::string CONSTANTOCEANBOUNDARY = "Nextsim::ConstantOceanBoundary";

template <>
Module<Nextsim::IOceanBoundary>::map Module<Nextsim::IOceanBoundary>::functionMap = {
    { CONSTANTOCEANBOUNDARY, newImpl<Nextsim::IOceanBoundary, Nextsim::ConstantOceanBoundary> },
};

template <>
Module<Nextsim::IOceanBoundary>::fn Module<Nextsim::IOceanBoundary>::spf
    = functionMap.at(CONSTANTOCEANBOUNDARY);
template <>
std::unique_ptr<Nextsim::IOceanBoundary> Module<Nextsim::IOceanBoundary>::staticInstance
    = std::move(newImpl<Nextsim::IOceanBoundary, Nextsim::ConstantOceanBoundary>());

template <> std::string Module<Nextsim::IOceanBoundary>::moduleName()
{
    return "Nextsim::IOceanBoundary";
}

template <> HelpMap& getHelpRecursive<Nextsim::IOceanBoundary>(HelpMap& map, bool getAll)
{
    const std::string& pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back(
        { pfx + "." + Module<Nextsim::IOceanBoundary>::moduleName(), ConfigType::MODULE,
            { CONSTANTOCEANBOUNDARY }, CONSTANTOCEANBOUNDARY, "", "MODULE DESCRIPTION HERE" });
    return map;
}
template <> Nextsim::IOceanBoundary& getImplementation<Nextsim::IOceanBoundary>()
{
    return getImplTemplate<Nextsim::IOceanBoundary, OceanBoundaryModule>();
}
template <> void setImplementation<Nextsim::IOceanBoundary>(const std::string& implName)
{
    setImplTemplate<OceanBoundaryModule>(implName);
}
template <> std::unique_ptr<Nextsim::IOceanBoundary> getInstance()
{
    return getInstTemplate<Nextsim::IOceanBoundary, OceanBoundaryModule>();
}
OceanBoundaryModule::Constructor OceanBoundaryModule::ctor;
OceanBoundaryModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IOceanBoundary, OceanBoundaryModule>();
}

} /* namespace Module */
