/*!
 * @file ModelComponent.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelComponent.hpp"

namespace Nextsim {

std::map<std::string, ModelComponent*> ModelComponent::registeredModules;
ModelArray* ModelComponent::sharedArrays[static_cast<size_t>(SharedArray::COUNT)];
const ModelArray* ModelComponent::protectedArrays[static_cast<size_t>(ProtectedArray::COUNT)];
std::map<ModelComponent::SharedArray, ModelArray*> ModelComponent::registeredArrays;
std::map<ModelComponent::SharedArray, std::set<ModelArray**>> ModelComponent::reservedArrays;
std::map<ModelComponent::SharedArray, std::set<const ModelArray**>>
    ModelComponent::reservedSemiArrays;
std::map<ModelComponent::ProtectedArray, const ModelArray*>
    ModelComponent::registeredProtectedArrays;
std::map<ModelComponent::ProtectedArray, std::set<const ModelArray**>>
    ModelComponent::reservedProtectedArrays;
ModelComponent::ModelComponent() { }

void ModelComponent::setAllModuleData(const ModelState& stateIn)
{
    for (auto entry : registeredModules) {
        entry.second->setData(stateIn);
    }
}
ModelState ModelComponent::getAllModuleState()
{
    ModelState overallState;
    for (auto entry : registeredModules) {
        overallState.merge(entry.second->getState());
    }
    return overallState;
}

void ModelComponent::registerModule() { registeredModules[getName()] = this; }

void ModelComponent::unregisterAllModules() { registeredModules.clear(); }

void ModelComponent::getAllFieldNames(
    std::set<std::string>& uF, std::set<std::string>& vF, std::set<std::string>& zF)
{
    for (auto entry : registeredModules) {
        uF.merge(entry.second->uFields());
        vF.merge(entry.second->vFields());
        zF.merge(entry.second->zFields());
    }
}

void ModelComponent::registerSharedArray(SharedArray type, ModelArray* addr)
{
    registeredArrays[type] = addr;
    for (ModelArray** addrAddr : reservedArrays[type]) {
        *addrAddr = addr;
    }
    for (const ModelArray** addrAddr : reservedSemiArrays[type]) {
        *addrAddr = addr;
    }

    // Assignment of pointer in array
    sharedArrays[static_cast<size_t>(type)] = addr;
}

void ModelComponent::requestSharedArray(SharedArray type, ModelArray** addr)
{
    if (registeredArrays.count(type) > 0) {
        *addr = registeredArrays[type];
    } else {
        reservedArrays[type].insert(addr);
    }
}

void ModelComponent::requestProtectedArray(SharedArray type, const ModelArray** addr)
{
    if (registeredArrays.count(type) > 0) {
        *addr = registeredArrays[type];
    } else {
        reservedSemiArrays[type].insert(addr);
    }
}

void ModelComponent::registerProtectedArray(ProtectedArray type, const ModelArray* addr)
{
    registeredProtectedArrays[type] = addr;
    for (const ModelArray** addrAddr : reservedProtectedArrays[type]) {
        *addrAddr = addr;
    }

    // Assignment of pointer in array
    protectedArrays[static_cast<size_t>(type)] = addr;
}

void ModelComponent::requestProtectedArray(ProtectedArray type, const ModelArray** addr)
{
    if (registeredProtectedArrays.count(type) > 0) {
        *addr = registeredProtectedArrays[type];
    } else {
        reservedProtectedArrays[type].insert(addr);
    }
}

} /* namespace Nextsim */
