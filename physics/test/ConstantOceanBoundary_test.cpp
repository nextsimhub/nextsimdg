/*!
 * @file ConstantOceanBoundary_test.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ConstantOceanBoundary.hpp"

#include "include/ModelArrayRef.hpp"
#include "include/ModelState.hpp"


namespace Nextsim {

TEST_SUITE_BEGIN("ConstantOceanBoundary");
TEST_CASE("ConstantOcean Qio calculation")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    ConstantOceanBoundary cob;

    cob.setData(ModelState::DataMap());
    cob.updateBefore(TimestepTime());
    ModelArrayRef<ModelComponent::SharedArray::Q_IO, MARBackingStore> qio(ModelComponent::getSharedArray());

    REQUIRE(qio[0] != 0.);
}
TEST_SUITE_END();
}
