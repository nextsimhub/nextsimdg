/*!
 * @file ERA5Atm_test.cpp
 *
 * @date Nov 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/TOPAZOcean.hpp"

#include "include/IFluxCalculation.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/Time.hpp"

#include <filesystem>

#define TO_STR(s) TO_STRI(s)
#define TO_STRI(s) #s
#ifndef TEST_FILE_SOURCE
#define TEST_FILE_SOURCE .
#endif

namespace Nextsim {

TEST_CASE("TOPAZOcean test", "[TOPAZOcean]")
{
    std::string filePath = "topaz_test128x128.nc";
    std::string sourceDir = TO_STR(TEST_FILE_SOURCE);
    // Copy the test file from the test source directory to the working directory
    if (!std::filesystem::exists(filePath)) {
        std::filesystem::copy(sourceDir + "/" + filePath, ".");
    }

    // In the real model, the array sizes will have been set by the restart file by this point
    size_t nx = 128;
    size_t ny = 128;
    size_t nxvertex = nx + 1;
    size_t nyvertex = ny + 1;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nxvertex);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, nyvertex);

    TOPAZOcean topaz;
    topaz.setFilePath(filePath);

    ModelArrayRef<ModelComponent::ProtectedArray::SST, MARConstBackingStore> sst(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::SSS, MARConstBackingStore> sss(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::MLD, MARConstBackingStore> mld(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::OCEAN_U, MARConstBackingStore> u(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::OCEAN_V, MARConstBackingStore> v(
        ModelComponent::getProtectedArray());

    TimePoint t1("2000-01-01T00:00:00Z");
    TimestepTime tst = { t1, Duration(600) };

    // Get the forcing fields at time 0
    topaz.updateBefore(tst);

    double mdi = -2.03703597633448608e90;

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032);
    REQUIRE(sst(45, 35) == -0.045035);
    REQUIRE(mld(45, 35) == 10.045035);

    TimePoint t2("2000-02-01T00:00:00Z");
    topaz.updateBefore({ t2, Duration(600) });

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032 - 1);
    REQUIRE(sst(45, 35) == -0.045035 - 1);
    REQUIRE(mld(45, 35) == 10.045035 + 1);

    TimePoint t12("2000-12-01T00:00:00Z");
    topaz.updateBefore({ t12, Duration(600) });

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032 - 11);
    REQUIRE(sst(45, 35) == -0.045035 - 11);
    REQUIRE(mld(45, 35) == 10.045035 + 11);

    // All times after the last time sample should use the last sample's data
    TimePoint t120("2010-01-01T00:00:00Z");
    topaz.updateBefore({ t120, Duration(600) });

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032 - 11);
    REQUIRE(sst(45, 35) == -0.045035 - 11);
    REQUIRE(mld(45, 35) == 10.045035 + 11);

    std::filesystem::remove(filePath);
}

}
