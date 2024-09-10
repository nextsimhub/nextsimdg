/*!
 * @file ModelComponent.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelComponent.hpp"

#include "include/MissingData.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

ModelComponent::ModelComponent()
{
    std::cout << "Constructing a ModelComponent" << std::endl;
}

void ModelComponent::setData(const ModelState::DataMap& state)
{
    if (state.count(maskName) > 0) {
        setOceanMask(state.at(maskName));
    }
}
/*
 * This assumes that the HField array size has already been set in the restart
 * reading routine. The mask, like all ModelArrays, is double precision,
 * where 0 (false) is land, >0 (true) is ocean.
 */
void ModelComponent::setOceanMask(const ModelArray& mask)
{
    // First check that the internal ocean mask has the correct allocated size. If not, then the
    // mask definitely needs to be assigned.
    bool maskMatch = mask.trueSize() == oceanMaskInternal().trueSize();

    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        if (!maskMatch)
            break;
        maskMatch = (mask[i] == oceanMaskInternal(i));
    }

    // If the old and new masks match, then return
    if (maskMatch)
        return;

    // Otherwise, set the new mask and (re)calculate the ocean point to grid point mapping.
    oceanMaskInternal() = mask;

    // Generate the oceanIndex to grid index mapping
    oceanIndex().clear();
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        if (oceanMaskInternal(i)) {
            oceanIndex().push_back(i);
        }
    }
}

// Fills the nOcean and OceanIndex variables for the zero land case
void ModelComponent::noLandMask()
{
    HField mask(ModelArray::Type::H);
    mask.resize();
    mask = 1;
    setOceanMask(mask);
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
    case (ModelArray::Type::Z): {
        ModelArray copy = ModelArray::ZField();
        copy = MissingData::value();
        size_t nZ = data.dimensions()[data.nDimensions() - 1];
        for (size_t i : oceanIndex()) {
            for (size_t k = 0; k < nZ; ++k) {
                copy.zIndexAndLayer(i, k) = data.zIndexAndLayer(i, k);
            }
        }
        return copy;
        break;
    }
    }
}

const ModelArray& ModelComponent::oceanMask() { return oceanMaskInternal(); }

} /* namespace Nextsim */
