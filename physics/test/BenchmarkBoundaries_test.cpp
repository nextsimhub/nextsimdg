/*!
 * @file BenchmarkBoundaries_test.cpp
 *
 * @date Sep 28, 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/BenchmarkAtmosphere.hpp"
#include "include/BenchmarkOcean.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("BenchmarkOcean");
TEST_CASE("OceanTest")
{
    // Expected dimensions of the benchmark domain
    const size_t nx = 256;
    const size_t ny = 256;
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

    double vMaxOcean = 0.01;
    REQUIRE(uOcean(0, 0) == -vMaxOcean);
    REQUIRE(vOcean(0, 0) == vMaxOcean);
    REQUIRE(uOcean(nx - 1, ny - 1) == (ny - 2.) / ny * vMaxOcean);
    REQUIRE(vOcean(nx - 1, ny - 1) == -(nx - 2.) / nx * vMaxOcean);
}

TEST_CASE("AtmosphereTest")
{
    // Expected dimensions of the benchmark domain
    const size_t nx = 256;
    const size_t ny = 256;
    ModelArray::setDimensions(ModelArray::Type::H, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::U, { nx, ny });
    ModelArray::setDimensions(ModelArray::Type::V, { nx, ny });

    BenchmarkAtmosphere benchAtm;
    benchAtm.setData(ModelState::DataMap());
    TimePoint time("2000-01-01T00:00:00");
    Duration step(3600);
    TimestepTime tst = { time, step };

    benchAtm.update(tst);

    // Get the u and v arrays
    ModelArrayRef<Protected::WIND_U> uWind(ModelComponent::getStore());
    ModelArrayRef<Protected::WIND_V> vWind(ModelComponent::getStore());
    // Check the wind at an arbitrary point lies in a reasonable range
    size_t iTest = 50;
    size_t jTest = 40;
    double uTest = uWind(iTest, jTest);
    double vTest = vWind(iTest, jTest);
    REQUIRE(uTest != 0.);
    REQUIRE(vTest != 0.);

    // Check that the cyclone is moving away (weakening) from the test point
    tst.start += tst.step;
    benchAtm.update(tst);
    REQUIRE(fabs(uWind(iTest, jTest) < fabs(uTest)));
    REQUIRE(fabs(vWind(iTest, jTest) < fabs(vTest)));
}

}
