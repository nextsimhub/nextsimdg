/*!
 * @file SlabOcean_test.cpp
 *
 * @date 10 Feb 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/SlabOcean.hpp"

#include "include/LinearFreezing.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("SlabOcean");
TEST_CASE("Test Qdw")
{
    // 1000 s timestep
    TimestepTime tst = { TimePoint(), Duration("P0-0T0:16:40") };
    double dt = tst.step.seconds();
    REQUIRE(dt == 1000);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    ModelArrayReferenceStore couplingArrays;

    double tOffset = 0.001;
    // Supply the data to the slab ocean
    HField sss(ModelArray::Type::H);
    sss = 32.;
    ModelComponent::getStore().registerArray(Protected::SSS, &sss, RO);

    HField sst(ModelArray::Type::H);
    sst = LinearFreezing()(sss[0]);
    ModelComponent::getStore().registerArray(Protected::SST, &sst, RO);

    HField mld(ModelArray::Type::H);
    mld = 6.48;
    ModelComponent::getStore().registerArray(Protected::MLD, &mld, RO);

    HField cpml(ModelArray::Type::H);
    cpml = Water::cp * Water::rho * mld;
    ModelComponent::getStore().registerArray(Protected::ML_BULK_CP, &cpml, RO);

    HField cice(ModelArray::Type::H);
    double cice0 = 0.5;
    cice = cice0;
    ModelComponent::getStore().registerArray(Protected::C_ICE, &cice, RO);

    HField data0(ModelArray::Type::H);
    data0 = 0;
    couplingArrays.registerArray(CouplingFields::Q_SS_NO_SW, &data0, RO);
    couplingArrays.registerArray(CouplingFields::Q_SS_SW, &data0, RO);
    couplingArrays.registerArray(CouplingFields::FWFLUX, &data0, RO);
    couplingArrays.registerArray(CouplingFields::SFLUX, &data0, RO);

    // External SS* data
    HField sssExt(ModelArray::Type::H);
    sssExt = sss;
    ModelComponent::getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);

    HField sstExt(ModelArray::Type::H);
    sstExt = sst + tOffset;
    ModelComponent::getStore().registerArray(Protected::EXT_SST, &sstExt, RO);

    SlabOcean slabOcean(couplingArrays);
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayRef<Protected::SLAB_QDW> qdw(ModelComponent::getStore());
    double prec = 1e-8;
    REQUIRE(qdw[0]
        == doctest::Approx(tOffset * cpml[0] / SlabOcean::defaultRelaxationTime).epsilon(prec));

    ModelArrayRef<Protected::SLAB_SST> sstSlab(ModelComponent::getStore());
    REQUIRE(sstSlab[0] != doctest::Approx(sst[0]).epsilon(prec / dt));
    REQUIRE(sstSlab[0] == doctest::Approx(sst[0] + dt * qdw[0] / cpml[0]).epsilon(prec));

    HField qswNet(ModelArray::Type::H);
    qswNet[0] = 15;
    couplingArrays.registerArray(CouplingFields::Q_SS_SW, &qswNet, RW);
    HField qNoSun(ModelArray::Type::H);
    qNoSun[0] = -17.5;
    couplingArrays.registerArray(CouplingFields::Q_SS_NO_SW, &qNoSun, RW);
    // Should not need to update anything else, as the slabOcean update only changes SLAB_SST
    slabOcean.update(tst);
    REQUIRE(sstSlab[0]
        == doctest::Approx(sst[0] - dt * (qswNet[0] + qNoSun[0] - qdw[0]) / cpml[0]).epsilon(prec));
}

TEST_CASE("Test Fdw")
{
    // 1000 s timestep
    TimestepTime tst = { TimePoint(), Duration("P0-0T0:16:40") };
    double dt = tst.step.seconds();
    REQUIRE(dt == 1000);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    ModelArrayReferenceStore couplingArrays;

    double sOffset = 0.1;
    // Supply the data to the slab ocean
    HField sss(ModelArray::Type::H);
    sss = 32.;
    ModelComponent::getStore().registerArray(Protected::SSS, &sss, RO);

    HField sst(ModelArray::Type::H);
    sst = LinearFreezing()(sss[0]);
    ModelComponent::getStore().registerArray(Protected::SST, &sst, RO);

    HField mld(ModelArray::Type::H);
    mld = 6.48;
    ModelComponent::getStore().registerArray(Protected::MLD, &mld, RO);

    HField cpml(ModelArray::Type::H);
    cpml = Water::cp * Water::rho * mld;
    ModelComponent::getStore().registerArray(Protected::ML_BULK_CP, &cpml, RO);

    HField data0(ModelArray::Type::H);
    data0 = 0;
    couplingArrays.registerArray(CouplingFields::Q_SS_NO_SW, &data0, RO);
    couplingArrays.registerArray(CouplingFields::Q_SS_SW, &data0, RO);
    couplingArrays.registerArray(CouplingFields::FWFLUX, &data0, RO);
    couplingArrays.registerArray(CouplingFields::SFLUX, &data0, RO);

    // External SS* data
    HField sssExt(ModelArray::Type::H);
    sssExt = sss + sOffset;
    ModelComponent::getStore().registerArray(Protected::EXT_SSS, &sssExt, RO);

    HField sstExt(ModelArray::Type::H);
    sstExt = sst;
    ModelComponent::getStore().registerArray(Protected::EXT_SST, &sstExt, RO);

    SlabOcean slabOcean(couplingArrays);
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayRef<Protected::SLAB_FDW> fdw(ModelComponent::getStore());
    double prec = 1e-6;
    REQUIRE(fdw[0]
        == doctest::Approx(
            -sOffset / sss[0] * mld[0] * Water::rho / SlabOcean::defaultRelaxationTime)
               .epsilon(prec));
    // Test that the finiteelement.cpp calculation of fdw is not being used
    double delS = -sOffset;
    double timeS = SlabOcean::defaultRelaxationTime;
    double ddt = tst.step.seconds();
    double oldFdw = delS * mld[0] * Water::rho / (timeS * sss[0] - ddt * delS);
    REQUIRE(fdw[0] != doctest::Approx(oldFdw).epsilon(prec * 1e-6));

    ModelArrayRef<Protected::SLAB_SSS> sssSlab(ModelComponent::getStore());
    REQUIRE(sssSlab[0] != doctest::Approx(sss[0]).epsilon(prec / dt));
    REQUIRE(sssSlab[0]
        == doctest::Approx(sss[0] - (fdw[0] * dt) / (mld[0] * Water::rho + fdw[0] * dt))
               .epsilon(prec));

    HField snowMeltFlux(ModelArray::Type::H);
    double snowMelt = -1e-4;
    double snowMeltVol = snowMelt * Ice::rhoSnow;
    snowMeltFlux = snowMeltVol / dt;
    couplingArrays.registerArray(CouplingFields::FWFLUX, &snowMeltFlux, RW);
    slabOcean.update(tst);
    REQUIRE(sssSlab[0]
        == doctest::Approx(sss[0]
            + (snowMeltVol - fdw[0] * dt) / (mld[0] * Water::rho - snowMeltVol + fdw[0] * dt))
               .epsilon(prec));
}
TEST_SUITE_END();
} /* namespace Nextsim */
