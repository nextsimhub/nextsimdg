/*!
 * @file IThermodynamicsModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IThermodynamicsModule.hpp"

#include "include/ThermoIce0.hpp"

#include <string>

namespace Module {
const std::string THERMOICE0 = "ThermoIce0";

template <>
Module<Nextsim::IThermodynamics>::map Module<Nextsim::IThermodynamics>::functionMap = {
    { THERMOICE0, newImpl<Nextsim::IThermodynamics, Nextsim::ThermoIce0> },
};

template <>
Module<Nextsim::IThermodynamics>::fn Module<Nextsim::IThermodynamics>::spf
    = functionMap.at(THERMOICE0);
template <>
std::unique_ptr<Nextsim::IThermodynamics> Module<Nextsim::IThermodynamics>::staticInstance
    = std::move(newImpl<Nextsim::IThermodynamics, Nextsim::ThermoIce0>());

template <> std::string Module<Nextsim::IThermodynamics>::moduleName() { return "IThermodynamics"; }

template <> Nextsim::IThermodynamics& getImplementation<Nextsim::IThermodynamics>()
{
    return getImplTemplate<Nextsim::IThermodynamics, IThermodynamicsModule>();
}
template <> void setImplementation<Nextsim::IThermodynamics>(const std::string& implName)
{
    setImplTemplate<IThermodynamicsModule>(implName);
}
template <> std::unique_ptr<Nextsim::IThermodynamics> getInstance()
{
    return getInstTemplate<Nextsim::IThermodynamics, IThermodynamicsModule>();
}
} /* namespace Module */
