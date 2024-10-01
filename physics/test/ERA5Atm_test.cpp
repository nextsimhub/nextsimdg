/*!
 * @file ERA5Atm_test.cpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "include/ERA5Atmosphere.hpp"

#include "include/IFluxCalculation.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"

#include <filesystem>
#include <memory>

#ifndef TEST_FILES_DIR
#define TEST_FILES_DIR "."
#endif

extern template class Module::Module<Nextsim::IFluxCalculation>;

namespace Nextsim {

// Null flux calculation
class NullFlux : public IFluxCalculation {
public:
    NullFlux()
    : IFluxCalculation()
    {
    }
    void update(const TimestepTime&) override { }

} nullFlux;

std::unique_ptr<IFluxCalculation> setNullFlux()
{
    return std::make_unique<NullFlux>();
}

TEST_SUITE_BEGIN("ERA5Atmosphere");
#ifdef USE_MPI
MPI_TEST_CASE("ERA5Atmosphere construction test", 1)
#else
TEST_CASE("ERA5Atmosphere construction test")
#endif
{
    const std::string filePath = "era5_test128x128.nc";
    const std::string orig_file = std::string(TEST_FILES_DIR) + "/" + filePath;
    std::filesystem::copy(orig_file, filePath, std::filesystem::copy_options::overwrite_existing);

    // In the real model, the array sizes will have been set by the restart file by this point
    size_t nx = 128;
    size_t ny = 128;
    size_t nxvertex = nx + 1;
    size_t nyvertex = ny + 1;

#ifdef USE_MPI
    ModelArray::setDimension(ModelArray::Dimension::X, nx, nx, 0);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nxvertex, nxvertex, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, nyvertex, nyvertex, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nxvertex);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, nyvertex);
#endif

    ERA5Atmosphere e5;

    Module::Module<IFluxCalculation>::setExternalImplementation(setNullFlux);

    e5.configure();
    e5.setFilePath(filePath);

    ModelArrayRef<Protected::T_AIR> tair(ModelComponent::getStore());
    ModelArrayRef<Protected::DEW_2M> tdew(ModelComponent::getStore());
    ModelArrayRef<Protected::P_AIR> pair(ModelComponent::getStore());
    ModelArrayRef<Protected::SW_IN> qswin(ModelComponent::getStore());
    ModelArrayRef<Protected::LW_IN> qlwin(ModelComponent::getStore());
    ModelArrayRef<Protected::WIND_SPEED> wind(ModelComponent::getStore());
    ModelArrayRef<Protected::WIND_U> u(ModelComponent::getStore());

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
    REQUIRE(u(30, 20) == (10 + 20.030) + 10 * 11);

    std::filesystem::remove(filePath);
}
TEST_SUITE_END();
}
