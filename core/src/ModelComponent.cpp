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

size_t ModelComponent::nOcean;
std::vector<size_t> ModelComponent::oceanIndex;

ModelComponent::ModelComponent() { noLandMask(); }

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
    nOcean = 0;
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

const ModelArray& ModelComponent::getOceanMask() { return oceanMaskSingleton(); }

} /* namespace Nextsim */
