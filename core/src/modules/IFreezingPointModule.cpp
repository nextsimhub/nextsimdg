/*!
 * @file IFreezingPointModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IFreezingPointModule.hpp"

#include "include/LinearFreezing.hpp"
#include "include/UnescoFreezing.hpp"

#include <string>

namespace Module {
const std::string LINEARFREEZING = "Nextsim::LinearFreezing";
const std::string UNESCOFREEZING = "Nextsim::UnescoFreezing";

template <>
Module<Nextsim::IFreezingPoint>::map Module<Nextsim::IFreezingPoint>::functionMap = {
    { LINEARFREEZING, newImpl<Nextsim::IFreezingPoint, Nextsim::LinearFreezing> },
    { UNESCOFREEZING, newImpl<Nextsim::IFreezingPoint, Nextsim::UnescoFreezing> },
};

template <>
Module<Nextsim::IFreezingPoint>::fn Module<Nextsim::IFreezingPoint>::spf
    = functionMap.at(LINEARFREEZING);
template <>
std::unique_ptr<Nextsim::IFreezingPoint> Module<Nextsim::IFreezingPoint>::staticInstance
    = std::move(newImpl<Nextsim::IFreezingPoint, Nextsim::LinearFreezing>());

template <> std::string Module<Nextsim::IFreezingPoint>::moduleName() { return "Nextsim::IFreezingPoint"; }

template <> HelpMap& getHelpRecursive<Nextsim::IFreezingPoint>(HelpMap& map, bool getAll)
{
    const std::string pfx = Nextsim::ConfiguredModule::MODULE_PREFIX;
    map[pfx].push_back({
        pfx + "." + Module<Nextsim::IFreezingPoint>::moduleName(), ConfigType::MODULE, {LINEARFREEZING, UNESCOFREEZING}, LINEARFREEZING, "",
                "The module selecting the model for the freezing point of sea water."
    });
    return map;
}
template <> Nextsim::IFreezingPoint& getImplementation<Nextsim::IFreezingPoint>()
{
    return getImplTemplate<Nextsim::IFreezingPoint, IFreezingPointModule>();
}
template <> void setImplementation<Nextsim::IFreezingPoint>(const std::string& implName)
{
    setImplTemplate<IFreezingPointModule>(implName);
}
template <> std::unique_ptr<Nextsim::IFreezingPoint> getInstance()
{
    return getInstTemplate<Nextsim::IFreezingPoint, IFreezingPointModule>();
}

IFreezingPointModule::Constructor IFreezingPointModule::ctor;
IFreezingPointModule::Constructor::Constructor()
{
   addToConfiguredModules<Nextsim::IFreezingPoint, IFreezingPointModule>();
}

} /* namespace Module */
