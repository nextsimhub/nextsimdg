/*!
 * @file IIceOceanHeatFluxModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceOceanHeatFluxModule.hpp"
#include "include/BasicIceOceanHeatFlux.hpp"
#include "include/ConfiguredIceOceanHeatFlux.hpp"

#include <string>

namespace Module {
const std::string BASICICEOCEANHEATFLUX = "Nextsim::BasicIceOceanHeatFlux";
const std::string CONFIGUREDICEOCEANHEATFLUX = "Nextsim::ConfiguredIceOceanHeatFlux";

template <>
Module<Nextsim::IIceOceanHeatFlux>::map Module<Nextsim::IIceOceanHeatFlux>::functionMap = {
    { BASICICEOCEANHEATFLUX, newImpl<Nextsim::IIceOceanHeatFlux, Nextsim::BasicIceOceanHeatFlux> },
    { CONFIGUREDICEOCEANHEATFLUX, newImpl<Nextsim::IIceOceanHeatFlux, Nextsim::ConfiguredIceOceanHeatFlux> },
};

template <>
Module<Nextsim::IIceOceanHeatFlux>::fn Module<Nextsim::IIceOceanHeatFlux>::spf
    = functionMap.at(BASICICEOCEANHEATFLUX);
template <>
std::unique_ptr<Nextsim::IIceOceanHeatFlux> Module<Nextsim::IIceOceanHeatFlux>::staticInstance
    = std::move(newImpl<Nextsim::IIceOceanHeatFlux, Nextsim::BasicIceOceanHeatFlux>());

template <> std::string Module<Nextsim::IIceOceanHeatFlux>::moduleName()
{
    return "Nextsim::IIceOceanHeatFlux";
}

template <> HelpMap& getHelpRecursive<Nextsim::IIceOceanHeatFlux>(HelpMap& map, bool getAll)
{
    const std::string pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({ pfx + "." + Module<Nextsim::IIceOceanHeatFlux>::moduleName(),
        ConfigType::MODULE, { BASICICEOCEANHEATFLUX, CONFIGUREDICEOCEANHEATFLUX }, BASICICEOCEANHEATFLUX, "",
        "The module for calculating the ice-ocean heat fluxes." });
    Nextsim::ConfiguredIceOceanHeatFlux::getHelpRecursive(map, getAll);
    return map;
}
template <> Nextsim::IIceOceanHeatFlux& getImplementation<Nextsim::IIceOceanHeatFlux>()
{
    return getImplTemplate<Nextsim::IIceOceanHeatFlux, IceOceanHeatFluxModule>();
}
template <> void setImplementation<Nextsim::IIceOceanHeatFlux>(const std::string& implName)
{
    setImplTemplate<IceOceanHeatFluxModule>(implName);
}
template <> std::unique_ptr<Nextsim::IIceOceanHeatFlux> getInstance()
{
    return getInstTemplate<Nextsim::IIceOceanHeatFlux, IceOceanHeatFluxModule>();
}
IceOceanHeatFluxModule::Constructor IceOceanHeatFluxModule::ctor;
IceOceanHeatFluxModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IIceOceanHeatFlux, IceOceanHeatFluxModule>();
}
} /* namespace Module */
