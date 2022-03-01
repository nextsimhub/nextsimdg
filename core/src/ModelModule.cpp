/*!
 * @file ModelModule.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelModule.hpp"

namespace Nextsim {

std::map<std::string, ModelModule*> ModelModule::registeredModules;

ModelModule::ModelModule()
{
    // TODO Auto-generated constructor stub
}

void ModelModule::setAllModuleData(const ModelState& stateIn)
{
    for (auto entry : registeredModules) {
        entry.second->setData(stateIn);
    }
}
ModelState ModelModule::getAllModuleState()
{
    ModelState overallState;
    for (auto entry : registeredModules) {
        overallState.merge(entry.second->getState());
    }
    return overallState;
}

void ModelModule::registerModule() { registeredModules[getName()] = this; }

void ModelModule::unregisterAllModules() { registeredModules.clear(); }
} /* namespace Nextsim */
