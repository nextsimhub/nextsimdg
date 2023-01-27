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

    ConstantOceanBoundary cob;

    cob.setData(ModelState::DataMap());

    ModelArrayRef<ModelComponent::SharedArray::Q_IO, MARBackingStore> qio(ModelComponent::getSharedArray());

    REQUIRE(qio[0] != 0.);
}
}
