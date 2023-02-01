/*!
 * @file SlabOcean_test.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/SlabOcean.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/LinearFreezing.hpp"
#include "include/constants.hpp"

namespace Nextsim {

// 1000 s timestep
TimestepTime tst = {TimePoint(0), Duration("P0-0T0:16:40")};

TEST_CASE("Test Qdw", "[SlabOcean]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    double tOffset = 0.001;

    // Supply the data to the slab ocean
    HField sss(ModelArray::Type::H);
    sss = 32.;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::SSS, &sss);

    HField sst(ModelArray::Type::H);
    sst = LinearFreezing()(sss[0]);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::SST, &sst);

    HField mld(ModelArray::Type::H);
    mld = 6.48;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::MLD, &mld);

    HField cpml(ModelArray::Type::H);
    cpml = Water::cp * Water::rho * mld;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::ML_BULK_CP, &cpml);

    HField data0(ModelArray::Type::H);
    data0 = 0;
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_IO, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_OW, &data0);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::C_ICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::DELTA_HICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::HSNOW_MELT, &data0);

    // External SS* data
    HField sssExt(ModelArray::Type::H);
    sssExt = sss;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EXT_SSS, &sssExt);

    HField sstExt(ModelArray::Type::H);
    sstExt = sst + tOffset;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EXT_SST, &sstExt);

    SlabOcean slabOcean;
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayRef<ModelComponent::ProtectedArray::SLAB_QDW, MARConstBackingStore> qdw(ModelComponent::getProtectedArray());
    REQUIRE(tst.step.seconds() == 1000);
    REQUIRE(qdw[0] == tOffset * cpml[0] * SlabOcean::defaultRelaxationTime);
}
} /* namespace Nextsim */
