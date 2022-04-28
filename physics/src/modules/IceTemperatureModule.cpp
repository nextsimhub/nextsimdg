/*!
 * @file IceTemperatureModule.cpp
 *
 * @date Apr 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IceTemperatureModule.hpp"

#include "include/ThermoIce0Temperature.hpp"

#include <string>

namespace Module {
const std::string CONSTANTICETEMPERATURE = "Nextsim::ConstantIceTemperature";
const std::string THERMOICE0TEMPERATURE = "Nextsim::ThermoIce0Temperature";

template <>
Module<Nextsim::IIceTemperature>::map Module<Nextsim::IIceTemperature>::functionMap = {
    { CONSTANTICETEMPERATURE, newImpl<Nextsim::IIceTemperature, Nextsim::ConstantIceTemperature> },
    { THERMOICE0TEMPERATURE, newImpl<Nextsim::IIceTemperature, Nextsim::ThermoIce0Temperature> },
};

template <>
Module<Nextsim::IIceTemperature>::fn Module<Nextsim::IIceTemperature>::spf
    = functionMap.at(CONSTANTICETEMPERATURE);
template <>
std::unique_ptr<Nextsim::IIceTemperature> Module<Nextsim::IIceTemperature>::staticInstance
    = std::move(newImpl<Nextsim::IIceTemperature, Nextsim::ConstantIceTemperature>());

template <> std::string Module<Nextsim::IIceTemperature>::moduleName()
{
    return "Nextsim::IIceTemperature";
}

template <> Nextsim::IIceTemperature& getImplementation<Nextsim::IIceTemperature>()
{
    return getImplTemplate<Nextsim::IIceTemperature, IceTemperatureModule>();
}
template <> void setImplementation<Nextsim::IIceTemperature>(const std::string& implName)
{
    setImplTemplate<IceTemperatureModule>(implName);
}
template <> std::unique_ptr<Nextsim::IIceTemperature> getInstance()
{
    return getInstTemplate<Nextsim::IIceTemperature, IceTemperatureModule>();
}
IceTemperatureModule::Constructor IceTemperatureModule::ctor;
IceTemperatureModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IIceTemperature, IceTemperatureModule>();
}

} /* namespace Module */
