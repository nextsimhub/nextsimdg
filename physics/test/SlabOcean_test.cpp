/*!
 * @file SlabOcean_test.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/SlabOcean.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/LinearFreezing.hpp"
#include "include/constants.hpp"

#include <iostream> //FIXME remove me

namespace Nextsim {

TEST_SUITE_BEGIN("SlabOcean");
TEST_CASE("Test Qdw")
{
    // 1000 s timestep
    TimestepTime tst = {TimePoint(), Duration("P0-0T0:16:40")};
    double dt = tst.step.seconds();
    REQUIRE(dt == 1000);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    double tOffset = 0.001;
    // Supply the data to the slab ocean
    HField sss(ModelArray::Type::H);
    sss = 32.;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::SSS, &sss);

    HField sst(ModelArray::Type::H);
    sst = LinearFreezing()(sss[0]);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::SST, &sst);

    HField mld(ModelArray::Type::H);
    mld = 6.48;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::MLD, &mld);

    HField cpml(ModelArray::Type::H);
    cpml = Water::cp * Water::rho * mld;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::ML_BULK_CP, &cpml);

    HField cice(ModelArray::Type::H);
    double cice0 = 0.5;
    cice = cice0;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::C_ICE, &cice);

    HField data0(ModelArray::Type::H);
    data0 = 0;
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_IO, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_OW, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::DELTA_HICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::NEW_ICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::HSNOW_MELT, &data0);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EVAP_MINUS_PRECIP, &data0);

    // External SS* data
    HField sssExt(ModelArray::Type::H);
    sssExt = sss;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EXT_SSS, &sssExt);

    HField sstExt(ModelArray::Type::H);
    sstExt = sst + tOffset;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EXT_SST, &sstExt);

    SlabOcean slabOcean;
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayRef<ModelComponent::ProtectedArray::SLAB_QDW, MARConstBackingStore> qdw(ModelComponent::getProtectedArray());
    double prec = 1e-6;
    REQUIRE(qdw[0] == doctest::Approx(tOffset * cpml[0] / SlabOcean::defaultRelaxationTime).epsilon(prec));

    ModelArrayRef<ModelComponent::ProtectedArray::SLAB_SST, MARConstBackingStore> sstSlab(ModelComponent::getProtectedArray());
    REQUIRE(sstSlab[0] != doctest::Approx(sst[0]).epsilon(prec / dt));
    REQUIRE(sstSlab[0] == doctest::Approx(sst[0] + dt * qdw[0] / cpml[0]).epsilon(prec));

    HField qow(ModelArray::Type::H);
    qow[0] = 15;
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_OW, &qow);
    HField qio(ModelArray::Type::H);
    qio[0] = -17.5;
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_IO, &qio);
    // Should not need to update anything else, as the slabOcean update only changes SLAB_SST
    slabOcean.update(tst);
    REQUIRE(sstSlab[0] == doctest::Approx(sst[0] - dt * (cice0 * qow[0] + cice0 * qio[0] - qdw[0]) / cpml[0]).epsilon(prec));
}

TEST_CASE("Test Fdw")
{
    // 1000 s timestep
    TimestepTime tst = {TimePoint(), Duration("P0-0T0:16:40")};
    double dt = tst.step.seconds();
    REQUIRE(dt == 1000);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    double sOffset = 0.1;
    // Supply the data to the slab ocean
    HField sss(ModelArray::Type::H);
    sss = 32.;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::SSS, &sss);

    HField sst(ModelArray::Type::H);
    sst = LinearFreezing()(sss[0]);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::SST, &sst);

    HField mld(ModelArray::Type::H);
    mld = 6.48;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::MLD, &mld);

    HField cpml(ModelArray::Type::H);
    cpml = Water::cp * Water::rho * mld;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::ML_BULK_CP, &cpml);

    HField data0(ModelArray::Type::H);
    data0 = 0;
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_IO, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_OW, &data0);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::C_ICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::DELTA_HICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::NEW_ICE, &data0);
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::HSNOW_MELT, &data0);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EVAP_MINUS_PRECIP, &data0);

    // External SS* data
    HField sssExt(ModelArray::Type::H);
    sssExt = sss + sOffset;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EXT_SSS, &sssExt);

    HField sstExt(ModelArray::Type::H);
    sstExt = sst;
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::EXT_SST, &sstExt);

    SlabOcean slabOcean;
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayRef<ModelComponent::ProtectedArray::SLAB_FDW, MARConstBackingStore> fdw(ModelComponent::getProtectedArray());
    double prec = 1e-6;
    REQUIRE(fdw[0] == doctest::Approx(-sOffset / sss[0] * mld[0] * Water::rho / SlabOcean::defaultRelaxationTime).epsilon(prec));
    // Test that the finiteelement.cpp calculation of fdw is not being used
    double delS = -sOffset;
    double timeS = SlabOcean::defaultRelaxationTime;
    double ddt = tst.step.seconds();
    double oldFdw = delS * mld[0] * Water::rho / (timeS * sss[0] - ddt * delS);
    REQUIRE(fdw[0] != doctest::Approx(oldFdw).epsilon(prec * 1e-6));

    ModelArrayRef<ModelComponent::ProtectedArray::SLAB_SSS, MARConstBackingStore> sssSlab(ModelComponent::getProtectedArray());
    REQUIRE(sssSlab[0] != doctest::Approx(sss[0]).epsilon(prec / dt));
    REQUIRE(sssSlab[0] == doctest::Approx(sss[0] - (fdw[0] * dt) / (mld[0] * Water::rho + fdw[0] * dt)).epsilon(prec));

    HField snowMelt(ModelArray::Type::H);
    snowMelt = -1e-4;
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::HSNOW_MELT, &snowMelt);
    slabOcean.update(tst);
    double snowMeltVol = snowMelt[0] * Ice::rhoSnow;
    REQUIRE(sssSlab[0] == doctest::Approx(sss[0] + (snowMeltVol - fdw[0] * dt) / (mld[0] * Water::rho - snowMeltVol + fdw[0] * dt)).epsilon(prec));
}
TEST_SUITE_END();
} /* namespace Nextsim */
