/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include <doctest/doctest.h>

#include "include/TOPAZOcean.hpp"

#include "include/IFluxCalculation.hpp"
#include "include/Time.hpp"

#include <filesystem>

#ifndef TEST_FILES_DIR
#define TEST_FILES_DIR "."
#endif

namespace Nextsim {

TEST_SUITE_BEGIN("TOPAZOcean");
TEST_CASE("TOPAZOcean test")
{
    const std::string filePath = "topaz_test128x128.nc";
    const std::string orig_file = std::string(TEST_FILES_DIR) + "/" + filePath;
    std::filesystem::copy(orig_file, filePath, std::filesystem::copy_options::overwrite_existing);

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
    topaz.configure();
    topaz.setFilePath(filePath);

    ModelArrayAccessor<Protected::EXT_SST> sstAccessor(ModelComponent::getStore());
    ModelArrayAccessor<Protected::EXT_SSS> sssAccessor(ModelComponent::getStore());
    ModelArrayAccessor<Protected::MLD> mldAccessor(ModelComponent::getStore());
    ModelArrayAccessor<Protected::OCEAN_U> uAccessor(ModelComponent::getStore());
    ModelArrayAccessor<Protected::OCEAN_V> vAccessor(ModelComponent::getStore());
    ModelArrayAccessor<Protected::SSH> sshAccessor(ModelComponent::getStore());

    TimePoint t1("2000-01-01T00:00:00Z");
    TimestepTime tst = { t1, Duration(600) };

    // The Qio calculation requires c_ice data
    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& cice = ciceAccessor.getHostRW();
    cice = 0.;

    // Get the forcing fields at time 0
    topaz.updateBefore(tst);

    FloatType mdi = -2.03703597633448608e90;

    // Use this, rather than the literal 0.035045, as the two are not equal at FloatType precision
    FloatType targetFrac = 35 * 0.001 + 45 * 0.000001;

    {
        const HField& sst = sstAccessor.getHostRO();
        const HField& sss = sssAccessor.getHostRO();
        const HField& mld = mldAccessor.getHostRO();
        const HField& u = uAccessor.getHostRO();
        const HField& v = vAccessor.getHostRO();
        const HField& ssh = sshAccessor.getHostRO();
        REQUIRE(sst(0, 0) == mdi);
        REQUIRE(sst(32, 32) == -0.032032_ft);
        REQUIRE(sst(45, 35) == -(0 + targetFrac));
        REQUIRE(mld(45, 35) == (10 + targetFrac));
        REQUIRE(ssh(45, 35) == (20 + targetFrac));
    }

    TimePoint t2("2000-02-01T00:00:00Z");
    topaz.updateBefore({ t2, Duration(600) });

    {
        const HField& sst = sstAccessor.getHostRO();
        const HField& sss = sssAccessor.getHostRO();
        const HField& mld = mldAccessor.getHostRO();
        const HField& u = uAccessor.getHostRO();
        const HField& v = vAccessor.getHostRO();
        const HField& ssh = sshAccessor.getHostRO();
        REQUIRE(sst(0, 0) == mdi);
        REQUIRE(sst(32, 32) == -0.032032_ft - 1);
        REQUIRE(sst(45, 35) == -(0 + targetFrac) - 1);
        REQUIRE(mld(45, 35) == (10 + targetFrac) + 1);
        REQUIRE(ssh(45, 35) == (20 + targetFrac) + 1);
    }

    TimePoint t12("2000-12-01T00:00:00Z");
    topaz.updateBefore({ t12, Duration(600) });

    {
        const HField& sst = sstAccessor.getHostRO();
        const HField& sss = sssAccessor.getHostRO();
        const HField& mld = mldAccessor.getHostRO();
        const HField& u = uAccessor.getHostRO();
        const HField& v = vAccessor.getHostRO();
        const HField& ssh = sshAccessor.getHostRO();
        REQUIRE(sst(0, 0) == mdi);
        REQUIRE(sst(32, 32) == -0.032032_ft - 11);
        REQUIRE(sst(45, 35) == -(0 + targetFrac) - 11);
        REQUIRE(mld(45, 35) == (10 + targetFrac) + 11);
        REQUIRE(ssh(45, 35) == (20 + targetFrac) + 11);
    }

    // All times after the last time sample should use the last sample's data
    TimePoint t120("2010-01-01T00:00:00Z");
    topaz.updateBefore({ t120, Duration(600) });

    {
        const HField& sst = sstAccessor.getHostRO();
        const HField& sss = sssAccessor.getHostRO();
        const HField& mld = mldAccessor.getHostRO();
        const HField& u = uAccessor.getHostRO();
        const HField& v = vAccessor.getHostRO();
        const HField& ssh = sshAccessor.getHostRO();
        REQUIRE(sst(0, 0) == mdi);
        REQUIRE(sst(32, 32) == -0.032032_ft - 11);
        REQUIRE(sst(45, 35) == -(0 + targetFrac) - 11);
        REQUIRE(mld(45, 35) == (10 + targetFrac) + 11);
        REQUIRE(ssh(45, 35) == (20 + targetFrac) + 11);
    }

    std::filesystem::remove(filePath);
}
TEST_SUITE_END();
}
