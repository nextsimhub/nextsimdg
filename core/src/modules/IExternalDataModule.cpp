/*!
 * @file IExternalDataModule.cpp
 *
 * @date Mar 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IExternalDataModule.hpp"

#include "include/ConstantExternalData.hpp"

#include <string>

namespace Module {
const std::string CONSTANTEXTERNALDATA = "Nextsim::ConstantExternalData";

template <>
Module<Nextsim::IExternalData>::map Module<Nextsim::IExternalData>::functionMap = {
    {CONSTANTEXTERNALDATA, newImpl<Nextsim::IExternalData, Nextsim::ConstantExternalData>},
};

template <>
Module<Nextsim::IExternalData>::fn Module<Nextsim::IExternalData>::spf = functionMap.at(CONSTANTEXTERNALDATA);
template <>
std::unique_ptr<Nextsim::IExternalData> Module<Nextsim::IExternalData>::staticInstance
= std::move(Module<Nextsim::IExternalData>::spf());

template <>
std::string Module<Nextsim::IExternalData>::moduleName(){    return "Nextsim::IExternalData";}

template <>
Nextsim::IExternalData& getImplementation<Nextsim::IExternalData>()
{
    return getImplTemplate<Nextsim::IExternalData, IExternalDataModule>();
}
template <>
void setImplementation<Nextsim::IExternalData>(const std::string& implName)
{
    setImplTemplate<IExternalDataModule>(implName);
}
template <>
std::unique_ptr<Nextsim::IExternalData> getInstance()
{
    return getInstTemplate<Nextsim::IExternalData, IExternalDataModule>();
}
IExternalDataModule::Constructor IExternalDataModule::ctor;
IExternalDataModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IExternalData, IExternalDataModule>();
}

} /* namespace Module */

} /* namespace Module */
