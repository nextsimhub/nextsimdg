/*!
 * @file BenchmarkBoundaries_test.cpp
 *
 * @date Sep 28, 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/BenchmarkOcean.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("BenchmarkOcean");
TEST_CASE("OceanTest")
{
    // Expected dimensions of the benchmark domain
    const size_t nx = 154;
    const size_t ny = 121;
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::U, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::V, { nx, ny });

    BenchmarkOcean benchOcean;
    benchOcean.setData(ModelState::DataMap());

    // Get the u and v arrays
    ModelArrayRef<Protected::OCEAN_U> uOcean(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_V> vOcean(ModelComponent::getStore());
    // Check the wind at an arbitrary point lies in a reasonable range
    size_t iTest = 50;
    size_t jTest = 40;
    REQUIRE(uOcean(iTest, jTest) != 0.);
    REQUIRE(vOcean(iTest, jTest) != 0.);

}

}
