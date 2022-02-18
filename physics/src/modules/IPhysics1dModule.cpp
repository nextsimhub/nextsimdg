/*!
 * @file IPhysics1dModule.cpp
 *
 * @date Feb 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IPhysics1dModule.hpp"

#include "include/NextsimPhysics.hpp"

namespace Module {

const std::string NEXTSIM_PHYSICS = "NextsimPhysics";

template <>
Module<Nextsim::IPhysics1d>::map Module<Nextsim::IPhysics1d>::functionMap = {
        {NEXTSIM_PHYSICS, newImpl<Nextsim::IPhysics1d, Nextsim::NextsimPhysics>},
};
template <>
Module<Nextsim::IPhysics1d>::fn Module<Nextsim::IPhysics1d>::spf = functionMap.at(NEXTSIM_PHYSICS);
template <>
std::unique_ptr<Nextsim::IPhysics1d> Module<Nextsim::IPhysics1d>::staticInstance = std::move(Module<Nextsim::IPhysics1d>::spf());

template <>
std::string Module<Nextsim::IPhysics1d>::moduleName()
{
    return "IPhysics1d";
}

template <>
Nextsim::IPhysics1d& getImplementation<Nextsim::IPhysics1d>()
{
    return getImplTemplate<Nextsim::IPhysics1d, IPhysics1dModule>();
}

template <>
void setImplementation<Nextsim::IPhysics1d>(const std::string& implName)
{
    setImplTemplate<IPhysics1dModule>(implName);
}

template <>
std::unique_ptr<Nextsim::IPhysics1d> getInstance()
{
    return getInstTemplate<Nextsim::IPhysics1d, IPhysics1dModule>();
}
} /* namespace Module */
