/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelComponent.hpp"

namespace Nextsim {

size_t ModelComponent::nOcean = 0;
std::vector<size_t> ModelComponent::oceanIndex;
bool ModelComponent::columnPhysicsStoreIsDestroyed; // initialized to 0 because of static lifetime

#if USE_KOKKOS
KokkosDeviceMapView<DeviceIndex> ModelComponent::oceanIndexDevice;
#endif

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

#if USE_KOKKOS
    makeOceanIndexDevice();
#endif
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

#if USE_KOKKOS
    makeOceanIndexDevice();
#endif
}

ModelArray ModelComponent::mask(const ModelArray& data, const double missingValue)
{
    auto copy = data;
    copy = missingValue;
    for (size_t iOcean = 0; iOcean < nOcean; ++iOcean) {
        copy(oceanIndex[iOcean]) = data(oceanIndex[iOcean]);
    }
    return copy;
}

const ModelArray& ModelComponent::oceanMask() { return oceanMaskSingleton(); }

#if USE_KOKKOS
void ModelComponent::makeOceanIndexDevice()
{
    // some tests don't need Kokkos
    if (!Kokkos::is_initialized()) {
        return;
    }

    std::vector<DeviceIndex> buf(oceanIndex.begin(), oceanIndex.end());
    oceanIndexDevice = makeKokkosDeviceViewMap("oceanIndex", buf, MakeViewOptions::ALWAYS_COPY);
    Finalizer::registerUnique(destroyOceanIndex);
}

void ModelComponent::destroyOceanIndex() { oceanIndexDevice.assign_data(nullptr); }
#endif

} /* namespace Nextsim */
