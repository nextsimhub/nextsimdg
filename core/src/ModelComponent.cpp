/*!
 * @file ModelComponent.cpp
 *
 * @date 02 Feb 2025
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
        return data * getOceanMask() + MissingData::value() * (1 - getOceanMask());
        break;
    }
    case (ModelArray::Type::Z): {
        ModelArray copy = ModelArray::ZField();
        copy = MissingData::value();
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

const ModelArray& ModelComponent::getOceanMask() { return *p_oceanMaskH; }

void ModelComponent::checkFields(const TimestepTime& tst)
{
    // numToCheck must be in scope for checkFieldElement as well
    for (numToCheck = 0; numToCheck < fieldsToCheck.size(); ++numToCheck) {

        const ModelArray* array = std::get<1>(fieldsToCheck[numToCheck]);

        int nLayers;
        if (array->nDimensions() == 3)
            nLayers = array->dimensions()[2];
        else if (array->nDimensions() == 2)
            nLayers = 1;
        else
            throw std::logic_error(
                "ModelComponent::checkFields expected a field with 2 or 3 dimensions.\n");

        // layerToCheck must be in scope for checkFieldElement as well
        for (layerToCheck = 0; layerToCheck < nLayers; ++layerToCheck)
            overElements(std::bind(&ModelComponent::checkFieldsElement, this, std::placeholders::_1,
                             std::placeholders::_2),
                tst);
    }
}

void ModelComponent::checkFieldsElement(size_t i, const TimestepTime& tst)
{
    const auto& x = fieldsToCheck[numToCheck];
    const double& value = std::get<1>(x)->zIndexAndLayer(i, layerToCheck);
    const auto& bounds = std::get<2>(x);
    if (value < bounds.first || value > bounds.second || std::isnan(value))
        throw std::runtime_error(std::get<0>(x) + " is out of bounds: " + std::to_string(value)
            + " is not within [" + std::to_string(bounds.first) + ","
            + std::to_string(bounds.second) + "].\n" + "Error at time step " + tst.start.format()
            + " and index " + std::to_string(i));
}

} /* namespace Nextsim */
