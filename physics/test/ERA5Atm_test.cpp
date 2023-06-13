/*!
 * @file ERA5Atm_test.cpp
 *
 * @date Nov 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ERA5Atmosphere.hpp"

#include "include/IFluxCalculation.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/Module.hpp"
#include "include/Time.hpp"

#include <filesystem>
#include <memory>

#define TO_STR(s) TO_STRI(s)
#define TO_STRI(s) #s
#ifndef TEST_FILE_SOURCE
#define TEST_FILE_SOURCE .
#endif

namespace Nextsim {

TEST_SUITE_BEGIN("ERA5Atmosphere");
TEST_CASE("ERA5Atmosphere construction test")
{
    std::string filePath = "era5_test128x128.nc";
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

    ERA5Atmosphere e5;

    // Null flux calculation
    class NullFlux : public IFluxCalculation {
    public:
        NullFlux()
            : IFluxCalculation()
        {
        }
        void update(const TimestepTime&) override { }

    } nullFlux;

    Module::Module<IFluxCalculation>::setExternalImplementation(
        []() { return std::make_unique<NullFlux>(); });

    e5.configure();
    e5.setFilePath(filePath);

    ModelArrayRef<ModelComponent::ProtectedArray::T_AIR, MARConstBackingStore> tair(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::DEW_2M, MARConstBackingStore> tdew(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::P_AIR, MARConstBackingStore> pair(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::SW_IN, MARConstBackingStore> qswin(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::LW_IN, MARConstBackingStore> qlwin(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::WIND_SPEED, MARConstBackingStore> wind(
        ModelComponent::getProtectedArray());
    ModelArrayRef<ModelComponent::ProtectedArray::WIND_U, MARConstBackingStore> u(
        ModelComponent::getProtectedArray());

    TimePoint t1("2000-01-01T00:00:00Z");
    TimestepTime tst = { t1, Duration(600) };


    // Get the forcing fields at time 0
    e5.update(tst);

    REQUIRE(wind(0, 0) == 0.);
    REQUIRE(wind(12, 12) == 12.012);
    REQUIRE(wind(30, 20) == 20.030);
    REQUIRE(pair(30, 20) == (1.01e5 + 20.030));

    TimePoint t2("2000-02-01T00:00:00Z");
    e5.update({ t2, Duration(600) });

    REQUIRE(wind(0, 0) == 0. + 100.);
    REQUIRE(wind(12, 12) == 12.012 + 100);
    REQUIRE(wind(30, 20) == 20.030 + 100);
    REQUIRE(pair(30, 20) == (1.01e5 + 20.030) + 1000);

    TimePoint t12("2000-12-01T00:00:00Z");
    e5.update({ t12, Duration(600) });

    REQUIRE(wind(0, 0) == 0. + 100. * 11);
    REQUIRE(wind(12, 12) == 12.012 + 100 * 11);
    REQUIRE(wind(30, 20) == 20.030 + 100 * 11);
    REQUIRE(pair(30, 20) == (1.01e5 + 20.030) + 1000 * 11);

    // All times after the last time sample should use the last sample's data
    TimePoint t120("2010-01-01T00:00:00Z");
    e5.update({ t120, Duration(600) });

    REQUIRE(wind(0, 0) == 0. + 100. * 11);
    REQUIRE(wind(12, 12) == 12.012 + 100 * 11);
    REQUIRE(wind(30, 20) == 20.030 + 100 * 11);
    REQUIRE(pair(30, 20) == (1.01e5 + 20.030) + 1000 * 11);
    REQUIRE(u(30, 20) == (10 + 30.020) + 10 * 11);

    std::filesystem::remove(filePath);
}
TEST_SUITE_END();
}
