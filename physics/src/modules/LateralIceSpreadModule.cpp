/*!
 * @file LateralIceSpreadModule.cpp
 *
 * @date Apr 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/LateralIceSpreadModule.hpp"

#include "include/HiblerSpread.hpp"

#include <string>

namespace Module {
const std::string HIBLERSPREAD = "Nextsim::HiblerSpread";

template <>
Module<Nextsim::ILateralIceSpread>::map Module<Nextsim::ILateralIceSpread>::functionMap = {
    { HIBLERSPREAD, newImpl<Nextsim::ILateralIceSpread, Nextsim::HiblerSpread> },
};

template <>
Module<Nextsim::ILateralIceSpread>::fn Module<Nextsim::ILateralIceSpread>::spf
    = functionMap.at(HIBLERSPREAD);
template <>
std::unique_ptr<Nextsim::ILateralIceSpread> Module<Nextsim::ILateralIceSpread>::staticInstance
    = std::move(newImpl<Nextsim::ILateralIceSpread, Nextsim::HiblerSpread>());

template <> std::string Module<Nextsim::ILateralIceSpread>::moduleName()
{
    return "Nextsim::ILateralIceSpread";
}

template <> Nextsim::ILateralIceSpread& getImplementation<Nextsim::ILateralIceSpread>()
{
    return getImplTemplate<Nextsim::ILateralIceSpread, LateralIceSpreadModule>();
}
template <> void setImplementation<Nextsim::ILateralIceSpread>(const std::string& implName)
{
    setImplTemplate<LateralIceSpreadModule>(implName);
}
template <> std::unique_ptr<Nextsim::ILateralIceSpread> getInstance()
{
    return getInstTemplate<Nextsim::ILateralIceSpread, LateralIceSpreadModule>();
}
LateralIceSpreadModule::Constructor LateralIceSpreadModule::ctor;
LateralIceSpreadModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::ILateralIceSpread, LateralIceSpreadModule>();
}

} /* namespace Module */
