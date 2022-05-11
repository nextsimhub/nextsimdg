/*!
 * @file OceanStateModule.cpp
 *
 * @date May 11, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/OceanStateModule.hpp"

#include "include/DummyOceanState.hpp"

#include <string>

namespace Module {
const std::string DUMMYOCEANSTATE = "Nextsim::DummyOceanState";

template <>
Module<Nextsim::OceanState>::map Module<Nextsim::OceanState>::functionMap = {
    {DUMMYOCEANSTATE, newImpl<Nextsim::OceanState, Nextsim::DummyOceanState>},
};

template <>
Module<Nextsim::OceanState>::fn Module<Nextsim::OceanState>::spf = functionMap.at(DUMMYOCEANSTATE);
template <>
std::unique_ptr<Nextsim::OceanState> Module<Nextsim::OceanState>::staticInstance
= std::move(newImpl<Nextsim::OceanState, Nextsim::DummyOceanState>());

template <>
std::string Module<Nextsim::OceanState>::moduleName(){    return "Nextsim::OceanState";}

template <>
Nextsim::OceanState& getImplementation<Nextsim::OceanState>()
{
    return getImplTemplate<Nextsim::OceanState, OceanStateModule>();
}
template <>
void setImplementation<Nextsim::OceanState>(const std::string& implName)
{
    setImplTemplate<OceanStateModule>(implName);
}
template <>
std::unique_ptr<Nextsim::OceanState> getInstance()
{
    return getInstTemplate<Nextsim::OceanState, OceanStateModule>();
}
OceanStateModule::Constructor OceanStateModule::ctor;
OceanStateModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::OceanState, OceanStateModule>();
}

} /* namespace Module */
