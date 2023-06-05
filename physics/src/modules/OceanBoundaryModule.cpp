/*!
 * @file OceanBoundaryModule.cpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/OceanBoundaryModule.hpp"

#include "include/ConstantOceanBoundary.hpp"
#include "include/ConfiguredOcean.hpp"
#include "include/FluxConfiguredOcean.hpp"
#include "include/IIceOceanHeatFlux.hpp"

#include <string>

namespace Module {
const std::string CONSTANTOCEANBOUNDARY = "Nextsim::ConstantOceanBoundary";
const std::string CONFIGUREDOCEAN = "Nextsim::ConfiguredOcean";
const std::string FLUXCONFIGUREDOCEAN = "Nextsim::FluxConfiguredOcean";

template <>
Module<Nextsim::IOceanBoundary>::map Module<Nextsim::IOceanBoundary>::functionMap = {
    { CONSTANTOCEANBOUNDARY, newImpl<Nextsim::IOceanBoundary, Nextsim::ConstantOceanBoundary> },
    { CONFIGUREDOCEAN, newImpl<Nextsim::IOceanBoundary, Nextsim::ConfiguredOcean> },
    { FLUXCONFIGUREDOCEAN, newImpl<Nextsim::IOceanBoundary, Nextsim::FluxConfiguredOcean> },
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
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IOceanBoundary>::moduleName(),
        ConfigType::MODULE, { CONSTANTOCEANBOUNDARY, CONFIGUREDOCEAN }, CONSTANTOCEANBOUNDARY, "",
        "Classes providing the oceanic inputs into the ice physics." });
    Nextsim::ConfiguredOcean::getHelpRecursive(map, getAll);
    Nextsim::FluxConfiguredOcean::getHelpRecursive(map, getAll);

    getHelpRecursive<Nextsim::IIceOceanHeatFlux>(map, getAll);
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
