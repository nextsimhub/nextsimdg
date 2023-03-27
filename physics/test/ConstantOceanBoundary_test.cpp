/*!
 * @file ConstantOceanBoundary_test.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/ConstantOceanBoundary.hpp"

#include "include/ModelArrayRef.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {

TEST_CASE("ConstantOcean Qio calculation", "[ConstantOceanBoundary]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    HField cice(ModelArray::Type::H);
    cice = 1.0; // Need some ice if Qio is to be calculated
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::C_ICE, &cice);
    ConstantOceanBoundary cob;

    cob.setData(ModelState::DataMap());
    cob.updateBefore(TimestepTime());
    ModelArrayRef<ModelComponent::SharedArray::Q_IO, MARBackingStore> qio(ModelComponent::getSharedArray());

    REQUIRE(qio[0] != 0.);
}
}
