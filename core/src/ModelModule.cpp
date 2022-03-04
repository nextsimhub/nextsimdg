/*!
 * @file ModelModule.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelModule.hpp"

namespace Nextsim {

std::map<std::string, ModelModule*> ModelModule::registeredModules;
std::map<ModelModule::SharedArray, ModelArray*> ModelModule::registeredArrays;
std::map<ModelModule::SharedArray, std::set<ModelArray**>> ModelModule::reservedArrays;
ModelModule::ModelModule() { }

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

void ModelModule::getAllFieldNames(
    std::set<std::string>& uF, std::set<std::string>& vF, std::set<std::string>& zF)
{
    for (auto entry : registeredModules) {
        uF.merge(entry.second->uFields());
        vF.merge(entry.second->vFields());
        zF.merge(entry.second->zFields());
    }
}

void ModelModule::registerSharedArray(SharedArray type, ModelArray* addr)
{
    registeredArrays[type] = addr;
    for (ModelArray** addrAddr : reservedArrays[type]) {
        *addrAddr = addr;
    }
}

void ModelModule::requestSharedArray(SharedArray type, ModelArray** addr)
{
    if (registeredArrays.count(type) > 0) {
        *addr = registeredArrays[type];
    } else {
        reservedArrays[type].insert(addr);
    }
}
} /* namespace Nextsim */
