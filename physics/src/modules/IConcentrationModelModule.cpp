/*!
 * @file IConcentrationModelModule.cpp
 *
 * @date Feb 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IConcentrationModelModule.hpp"

#include "include/HiblerConcentration.hpp"

#include <string>

namespace Module {
const std::string HIBLERCONCENTRATION = "HiblerConcentration";

template <>
Module<Nextsim::IConcentrationModel>::map Module<Nextsim::IConcentrationModel>::functionMap = {
    { HIBLERCONCENTRATION, newImpl<Nextsim::IConcentrationModel, Nextsim::HiblerConcentration> },
};

template <>
Module<Nextsim::IConcentrationModel>::fn Module<Nextsim::IConcentrationModel>::spf
    = functionMap.at(HIBLERCONCENTRATION);
template <>
std::unique_ptr<Nextsim::IConcentrationModel> Module<Nextsim::IConcentrationModel>::staticInstance
    = std::move(newImpl<Nextsim::IConcentrationModel, Nextsim::HiblerConcentration>());

template <> std::string Module<Nextsim::IConcentrationModel>::moduleName()
{
    return "IConcentrationModel";
}

template <> Nextsim::IConcentrationModel& getImplementation<Nextsim::IConcentrationModel>()
{
    return getImplTemplate<Nextsim::IConcentrationModel, IConcentrationModelModule>();
}
template <> void setImplementation<Nextsim::IConcentrationModel>(const std::string& implName)
{
    setImplTemplate<IConcentrationModelModule>(implName);
}
template <> std::unique_ptr<Nextsim::IConcentrationModel> getInstance()
{
    return getInstTemplate<Nextsim::IConcentrationModel, IConcentrationModelModule>();
}
} /* namespace Module */
