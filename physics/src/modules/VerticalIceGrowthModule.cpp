/*!
 * @file VerticalIceGrowthModule.cpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/VerticalIceGrowthModule.hpp"

#include "include/ThermoIce0Growth.hpp"

#include <string>

namespace Module {
const std::string THERMOICE0GROWTH = "Nextsim::ThermoIce0Growth";

template <>
Module<Nextsim::IVerticalIceGrowth>::map Module<Nextsim::IVerticalIceGrowth>::functionMap = {
    { THERMOICE0GROWTH, newImpl<Nextsim::IVerticalIceGrowth, Nextsim::ThermoIce0Growth> },
};

template <>
Module<Nextsim::IVerticalIceGrowth>::fn Module<Nextsim::IVerticalIceGrowth>::spf
    = functionMap.at(THERMOICE0GROWTH);
template <>
std::unique_ptr<Nextsim::IVerticalIceGrowth> Module<Nextsim::IVerticalIceGrowth>::staticInstance
    = std::move(Module<Nextsim::IVerticalIceGrowth>::spf());

template <> std::string Module<Nextsim::IVerticalIceGrowth>::moduleName()
{
    return "Nextsim::IVerticalIceGrowth";
}

template <> Nextsim::IVerticalIceGrowth& getImplementation<Nextsim::IVerticalIceGrowth>()
{
    return getImplTemplate<Nextsim::IVerticalIceGrowth, VerticalIceGrowthModule>();
}
template <> void setImplementation<Nextsim::IVerticalIceGrowth>(const std::string& implName)
{
    setImplTemplate<VerticalIceGrowthModule>(implName);
}
template <> std::unique_ptr<Nextsim::IVerticalIceGrowth> getInstance()
{
    return getInstTemplate<Nextsim::IVerticalIceGrowth, VerticalIceGrowthModule>();
}
VerticalIceGrowthModule::Constructor VerticalIceGrowthModule::ctor;
VerticalIceGrowthModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IVerticalIceGrowth, VerticalIceGrowthModule>();
}

} /* namespace Module */
