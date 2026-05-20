/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
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
    ModelArrayAccessor<Protected::OCEAN_U> uOceanAccessor(ModelComponent::getStore());
    const UField& uOcean = uOceanAccessor.getHostRO();
    ModelArrayAccessor<Protected::OCEAN_V> vOceanAccessor(ModelComponent::getStore());
    const VField& vOcean = vOceanAccessor.getHostRO();
    // Check the wind at an arbitrary point lies in a reasonable range
    size_t iTest = 50;
    size_t jTest = 40;
    REQUIRE(uOcean(iTest, jTest) != 0.);
    REQUIRE(vOcean(iTest, jTest) != 0.);

    FloatType vMaxOcean = 0.01;
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
    ModelArrayAccessor<Protected::WIND_U> uWindAccessor(ModelComponent::getStore());
    ModelArrayAccessor<Protected::WIND_V> vWindAccessor(ModelComponent::getStore());
    // Check the wind at an arbitrary point lies in a reasonable range
    size_t iTest = 50;
    size_t jTest = 40;
    const FloatType uTest = uWindAccessor.getHostRO()(iTest, jTest);
    const FloatType vTest = vWindAccessor.getHostRO()(iTest, jTest);
    REQUIRE(uTest != 0.);
    REQUIRE(vTest != 0.);

    // Check that the cyclone is moving away (weakening) from the test point
    tst.start += tst.step;
    benchAtm.update(tst);
    REQUIRE(fabs(uWindAccessor.getHostRO()(iTest, jTest) < fabs(uTest)));
    REQUIRE(fabs(vWindAccessor.getHostRO()(iTest, jTest) < fabs(vTest)));
}

}
