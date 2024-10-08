/*!
 * @file ERA5Atm_test.cpp
 *
 * @date 23 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "include/TOPAZOcean.hpp"

#include "include/IFluxCalculation.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/Time.hpp"

#include <filesystem>

#ifndef TEST_FILES_DIR
#define TEST_FILES_DIR "."
#endif

namespace Nextsim {

TEST_SUITE_BEGIN("TOPAZOcean");
#ifdef USE_MPI
MPI_TEST_CASE("TOPAZOcean test", 1)
#else
TEST_CASE("TOPAZOcean test")
#endif
{
    const std::string filePath = "topaz_test128x128.nc";
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

    TOPAZOcean topaz;
    topaz.configure();
    topaz.setFilePath(filePath);

    ModelArrayRef<Protected::EXT_SST> sst(ModelComponent::getStore());
    ModelArrayRef<Protected::EXT_SSS> sss(ModelComponent::getStore());
    ModelArrayRef<Protected::MLD> mld(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_U> u(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_V> v(ModelComponent::getStore());
    ModelArrayRef<Protected::SSH> ssh(ModelComponent::getStore());

    TimePoint t1("2000-01-01T00:00:00Z");
    TimestepTime tst = { t1, Duration(600) };

    // The Qio calculation requires c_ice data
    HField cice(ModelArray::Type::H);
    cice = 0.;
    ModelComponent::getStore().registerArray(Protected::C_ICE, &cice, RO);

    // Get the forcing fields at time 0
    topaz.updateBefore(tst);

    double mdi = -2.03703597633448608e90;

    // Use this, rather than the literal 0.035045, as the two are not equal at double precision
    double targetFrac = 35 * 0.001 + 45 * 0.000001;

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032);
    REQUIRE(sst(45, 35) == -(0 + targetFrac));
    REQUIRE(mld(45, 35) == (10 + targetFrac));
    REQUIRE(ssh(45, 35) == (20 + targetFrac));

    TimePoint t2("2000-02-01T00:00:00Z");
    topaz.updateBefore({ t2, Duration(600) });

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032 - 1);
    REQUIRE(sst(45, 35) == -(0 + targetFrac) - 1);
    REQUIRE(mld(45, 35) == (10 + targetFrac) + 1);
    REQUIRE(ssh(45, 35) == (20 + targetFrac) + 1);

    TimePoint t12("2000-12-01T00:00:00Z");
    topaz.updateBefore({ t12, Duration(600) });

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032 - 11);
    REQUIRE(sst(45, 35) == -(0 + targetFrac) - 11);
    REQUIRE(mld(45, 35) == (10 + targetFrac) + 11);
    REQUIRE(ssh(45, 35) == (20 + targetFrac) + 11);

    // All times after the last time sample should use the last sample's data
    TimePoint t120("2010-01-01T00:00:00Z");
    topaz.updateBefore({ t120, Duration(600) });

    REQUIRE(sst(0, 0) == mdi);
    REQUIRE(sst(32, 32) == -0.032032 - 11);
    REQUIRE(sst(45, 35) == -(0 + targetFrac) - 11);
    REQUIRE(mld(45, 35) == (10 + targetFrac) + 11);
    REQUIRE(ssh(45, 35) == (20 + targetFrac) + 11);

    std::filesystem::remove(filePath);
}
TEST_SUITE_END();
}
