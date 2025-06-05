/*!
 * @file ModelComponent.cpp
 *
 * @date 25 Apr 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelComponent.hpp"

#include "include/MissingData.hpp"

namespace Nextsim {

size_t ModelComponent::nOcean = 0;
std::vector<size_t> ModelComponent::oceanIndex;

ModelComponent::ModelComponent()
{
    // We only set no land mask if the mask hasn't been set by someone else.
    if (nOcean == 0)
        noLandMask();
}

/*
 * This assumes that the HField array size has already been set in the restart
 * reading routine. The mask, like all ModelArrays, is double precision,
 * where 0 (false) is land, >0 (true) is ocean.
 */
void ModelComponent::setOceanMask(const ModelArray& mask)
{
    oceanMaskSingleton() = mask;
    // Generate the oceanIndex to grid index mapping
    // 1. Count the number of non-land squares
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        if (oceanMask()[i] > 0)
            ++nOcean;
    }
    oceanIndex.resize(nOcean);
    size_t iOceanIndex = 0;
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        if (oceanMask()[i] > 0) {
            oceanIndex[iOceanIndex++] = i;
        }
    }
}

// Fills the nOcean and OceanIndex variables for the zero land case
void ModelComponent::noLandMask()
{
    ModelArray newOceanMask(ModelArray::Type::H);
    newOceanMask.resize();
    newOceanMask = 1.; // All ocean
    oceanMaskSingleton() = newOceanMask;

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
        return data * oceanMask() + MissingData::value() * (1 - oceanMask());
        break;
    }
    }
}

const ModelArray& ModelComponent::oceanMask() { return oceanMaskSingleton(); }

} /* namespace Nextsim */
