/*!
 * @file ModelComponent.cpp
 *
 * @date 28 Feb 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelComponent.hpp"

#include "include/MissingData.hpp"

#include <cmath>

namespace Nextsim {

std::unordered_map<std::string, ModelComponent*> ModelComponent::registeredModules;
ModelArrayReferenceStore ModelComponent::store;
ModelArray* ModelComponent::p_oceanMaskH = nullptr;
size_t ModelComponent::nOcean;
std::vector<size_t> ModelComponent::oceanIndex;

ModelComponent::ModelComponent() { noLandMask(); }

void ModelComponent::setAllModuleData(const ModelState& stateIn)
{
    for (auto entry : registeredModules) {
        entry.second->setData(stateIn.data);
    }
}
ModelState ModelComponent::getAllModuleState()
{
    ModelState overallState;
    for (auto entry : registeredModules) {
        overallState.data.merge(entry.second->getState().data);
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

/*
 * This assumes that the HField array size has already been set in the restart
 * reading routine. The mask, like all ModelArrays, is double precision,
 * where 0 (false) is land, >0 (true) is ocean.
 */
void ModelComponent::setOceanMask(const ModelArray& mask)
{
    if (p_oceanMaskH)
        delete p_oceanMaskH;
    p_oceanMaskH = new ModelArray(ModelArray::Type::H);
    ModelArray& oceanMaskH = *p_oceanMaskH;
    oceanMaskH.resize();
    oceanMaskH = mask;

    // Generate the oceanIndex to grid index mapping
    // 1. Count the number of non-land squares
    nOcean = 0;
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
    if (p_oceanMaskH)
        delete p_oceanMaskH;
    p_oceanMaskH = new ModelArray(ModelArray::Type::H);
    p_oceanMaskH->resize();
    *p_oceanMaskH = 1.; // All ocean

    nOcean = ModelArray::size(ModelArray::Type::H);
    oceanIndex.resize(nOcean);
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        oceanIndex[i] = i;
    }
}

ModelArray ModelComponent::mask(const ModelArray& data, const double missingValue)
{
    auto copy = data;
    copy = missingValue;
    const size_t nZ = data.getType() == ModelArray::Type::Z ? data.dimensions()[2] : 1;
    for (size_t iOcean = 0; iOcean < nOcean; ++iOcean) {
        const size_t i = oceanIndex[iOcean];
        for (size_t k = 0; k < nZ; ++k) {
            copy.zIndexAndLayer(i, k) = data.zIndexAndLayer(i, k);
        }
    }
    return copy;
}

const ModelArray& ModelComponent::getOceanMask() { return *p_oceanMaskH; }

} /* namespace Nextsim */
