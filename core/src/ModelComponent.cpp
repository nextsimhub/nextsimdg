/*!
 * @file ModelComponent.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelComponent.hpp"

namespace Nextsim {

std::unordered_map<std::string, ModelComponent*> ModelComponent::registeredModules;
ModelArray* ModelComponent::sharedArrays[static_cast<size_t>(SharedArray::COUNT)];
const ModelArray* ModelComponent::protectedArrays[static_cast<size_t>(ProtectedArray::COUNT)];
ModelArray ModelComponent::oceanMaskH = ModelArray::HField();
size_t ModelComponent::nOcean;
std::vector<size_t> ModelComponent::oceanIndex;

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

void ModelComponent::getAllFieldNames(std::unordered_set<std::string>& uF,
    std::unordered_set<std::string>& vF, std::unordered_set<std::string>& zF)
{
    for (auto entry : registeredModules) {
        uF.merge(entry.second->uFields());
        vF.merge(entry.second->vFields());
        zF.merge(entry.second->zFields());
    }
}

void ModelComponent::registerSharedArray(SharedArray type, ModelArray* addr)
{
    // Assignment of pointer in array
    sharedArrays[static_cast<size_t>(type)] = addr;
}

void ModelComponent::registerProtectedArray(ProtectedArray type, const ModelArray* addr)
{
    // Assignment of pointer in array
    protectedArrays[static_cast<size_t>(type)] = addr;
}

/*
 * This assumes that the HField array size has already been set in the restart
 * reading routine. The mask, like all ModelArrays, is double precision,
 * where 0 (false) is land, >0 (true) is ocean.
 */
void ModelComponent::setOceanMask(const ModelArray& mask)
{
    oceanMaskH.resize();
    oceanMaskH = mask;

    // Generate the oceanIndex to grid index mapping
    // 1. Count the number of non-land squares
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        if (oceanMaskH[i] > 0)
            ++nOcean;
    }
    oceanIndex.resize(nOcean);
    size_t iOceanIndex = 0;
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        if (oceanMaskH[i] > 0) {
            oceanIndex[iOceanIndex++] = i;
        }
    }
}

// Fills the nOcean and OceanIndex variables for the zero land case
void ModelComponent::noLandMask()
{
    nOcean = ModelArray::size(ModelArray::Type::H);
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        oceanIndex[i] = i;
    }
}

ModelArray ModelComponent::mask(const ModelArray& data)
{
    switch (data.getType()) {
    default: {
        return ModelArray(data);
        break;
    }
    case (ModelArray::Type::H):
    case (ModelArray::Type::U):
    case (ModelArray::Type::V): {
        return data * (oceanMaskH) + MissingData::value * (1 - oceanMaskH);
        break;
    }
    case (ModelArray::Type::Z): {
        ModelArray copy = ModelArray::ZField();
        copy = MissingData::value;
        size_t nZ = data.dimensions()[data.nDimensions() - 1];
        for (size_t iOcean = 0; iOcean < nOcean; ++iOcean) {
            size_t i = oceanIndex[iOcean];
            for (size_t k = 0; k < nZ; ++k) {
                copy.zIndexAndLayer(i, k) = data.zIndexAndLayer(i, k);
            }
        }
        return copy;
        break;
    }
    }
}

} /* namespace Nextsim */
