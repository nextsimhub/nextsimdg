/*!
 * @file AtmosphereBoundaryModule.cpp
 *
 * @date Sep 23, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/AtmosphereBoundaryModule.hpp"

#include "include/ConfiguredAtmosphere.hpp"
#include "include/ConstantAtmosphereBoundary.hpp"
#include "include/FluxConfiguredAtmosphere.hpp"
#include "include/MU71Atmosphere.hpp"

#include <string>

namespace Module {
const std::string CONSTANTATMOSPHEREBOUNDARY = "Nextsim::ConstantAtmosphereBoundary";
const std::string CONFIGUREDATMOSPHERE = "Nextsim::ConfiguredAtmosphere";
const std::string FLUXCONFIGUREDATMOSPHERE = "Nextsim::FluxConfiguredAtmosphere";
const std::string MU71ATMOSPHERE = "Nextsim::MU71Atmosphere";

template <>
Module<Nextsim::IAtmosphereBoundary>::map Module<Nextsim::IAtmosphereBoundary>::functionMap = {
    { CONSTANTATMOSPHEREBOUNDARY,
        newImpl<Nextsim::IAtmosphereBoundary, Nextsim::ConstantAtmosphereBoundary> },
    { CONFIGUREDATMOSPHERE, newImpl<Nextsim::IAtmosphereBoundary, Nextsim::ConfiguredAtmosphere> },
    { FLUXCONFIGUREDATMOSPHERE,
        newImpl<Nextsim::IAtmosphereBoundary, Nextsim::FluxConfiguredAtmosphere> },
    { MU71ATMOSPHERE, newImpl<Nextsim::IAtmosphereBoundary, Nextsim::MU71Atmosphere> },
};

template <>
Module<Nextsim::IAtmosphereBoundary>::fn Module<Nextsim::IAtmosphereBoundary>::spf
    = functionMap.at(CONSTANTATMOSPHEREBOUNDARY);
template <>
std::unique_ptr<Nextsim::IAtmosphereBoundary> Module<Nextsim::IAtmosphereBoundary>::staticInstance
    = std::move(newImpl<Nextsim::IAtmosphereBoundary, Nextsim::ConstantAtmosphereBoundary>());

template <> std::string Module<Nextsim::IAtmosphereBoundary>::moduleName()
{
    return "Nextsim::IAtmosphereBoundary";
}

template <> HelpMap& getHelpRecursive<Nextsim::IAtmosphereBoundary>(HelpMap& map, bool getAll)
{
    const std::string& pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IAtmosphereBoundary>::moduleName(),
        ConfigType::MODULE, { CONSTANTATMOSPHEREBOUNDARY }, CONSTANTATMOSPHEREBOUNDARY, "",
        "A Module to provide atmospheric inputs to the model." });
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IAtmosphereBoundary>::moduleName(),
        ConfigType::MODULE, { CONFIGUREDATMOSPHERE }, CONFIGUREDATMOSPHERE, "",
        "A configurable Module to provide atmospheric inputs to the model based on prescribed "
        "atmospheric state." });
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IAtmosphereBoundary>::moduleName(),
        ConfigType::MODULE, { FLUXCONFIGUREDATMOSPHERE }, FLUXCONFIGUREDATMOSPHERE, "",
        "A configurable Module to provide atmospheric inputs to the model based on prescribed "
        "atmospheric fluxes." });
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IAtmosphereBoundary>::moduleName(),
        ConfigType::MODULE, { MU71ATMOSPHERE }, MU71ATMOSPHERE, "",
        "A Module to provide atmospheric inputs to the model based on prescribed atmospheric "
        "fluxes according to Maykut and Untersteiner (1971)." });
    Nextsim::ConfiguredAtmosphere::getHelpRecursive(map, getAll);
    Nextsim::FluxConfiguredAtmosphere::getHelpRecursive(map, getAll);
    return map;
}
template <> Nextsim::IAtmosphereBoundary& getImplementation<Nextsim::IAtmosphereBoundary>()
{
    return getImplTemplate<Nextsim::IAtmosphereBoundary, AtmosphereBoundaryModule>();
}
template <> void setImplementation<Nextsim::IAtmosphereBoundary>(const std::string& implName)
{
    setImplTemplate<AtmosphereBoundaryModule>(implName);
}
template <> std::unique_ptr<Nextsim::IAtmosphereBoundary> getInstance()
{
    return getInstTemplate<Nextsim::IAtmosphereBoundary, AtmosphereBoundaryModule>();
}
AtmosphereBoundaryModule::Constructor AtmosphereBoundaryModule::ctor;
AtmosphereBoundaryModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IAtmosphereBoundary, AtmosphereBoundaryModule>();
}

} /* namespace Module */
