/*!
 * @file IceThermodynamicsModule.cpp
 *
 * @date Mar 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceThermodynamicsModule.hpp"

#include "include/ThermoIce0.hpp"
#include <string>

namespace Module {
const std::string THERMOICE0GROWTH = "Nextsim::ThermoIce0Growth";

template <>
Module<Nextsim::IIceThermodynamics>::map Module<Nextsim::IIceThermodynamics>::functionMap = {
    { THERMOICE0GROWTH, newImpl<Nextsim::IIceThermodynamics, Nextsim::ThermoIce0> },
};

template <>
Module<Nextsim::IIceThermodynamics>::fn Module<Nextsim::IIceThermodynamics>::spf
    = functionMap.at(THERMOICE0GROWTH);
template <>
std::unique_ptr<Nextsim::IIceThermodynamics> Module<Nextsim::IIceThermodynamics>::staticInstance
    = std::move(newImpl<Nextsim::IIceThermodynamics, Nextsim::ThermoIce0>());

template <> std::string Module<Nextsim::IIceThermodynamics>::moduleName()
{
    return "Nextsim::IVerticalIceGrowth";
}

template <> Nextsim::IIceThermodynamics& getImplementation<Nextsim::IIceThermodynamics>()
{
    return getImplTemplate<Nextsim::IIceThermodynamics, IceThermodynamicsModule>();
}
template <> void setImplementation<Nextsim::IIceThermodynamics>(const std::string& implName)
{
    setImplTemplate<IceThermodynamicsModule>(implName);
}
template <> std::unique_ptr<Nextsim::IIceThermodynamics> getInstance()
{
    return getInstTemplate<Nextsim::IIceThermodynamics, IceThermodynamicsModule>();
}
IceThermodynamicsModule::Constructor IceThermodynamicsModule::ctor;
IceThermodynamicsModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IIceThermodynamics, IceThermodynamicsModule>();
}

} /* namespace Module */
