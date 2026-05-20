/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include <doctest/doctest.h>

#include "include/SlabOcean.hpp"

#include "include/LinearFreezing.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("SlabOcean");
TEST_CASE("Test Qdw")
{
    // 1000 s timestep
    TimestepTime tst = { TimePoint(), Duration("P0-0T0:16:40") };
    FloatType dt = tst.step.seconds();
    REQUIRE(dt == 1000);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    ModelArrayStore couplingArrays;

    FloatType tOffset = 0.001;
    // Supply the data to the slab ocean
    ModelArrayAccessor<Protected::SSS, RW> sssAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sss = sssAccessor.getHostRW();
    sss = 32.;

    ModelArrayAccessor<Protected::SST, RW> sstAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sst = sstAccessor.getHostRW();
    sst = LinearFreezing()(sss[0]);

    ModelArrayAccessor<Protected::MLD, RW> mldAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& mld = mldAccessor.getHostRW();
    mld = 6.48;

    ModelArrayAccessor<Protected::ML_BULK_CP, RW> cpmlAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& cpml = cpmlAccessor.getHostRW();
    cpml = Water::cp * Water::rho * mld;

    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& cice = ciceAccessor.getHostRW();
    FloatType cice0 = 0.5;
    cice = cice0;

    /*
    ModelArrayAccessor<CouplingFields::Q_SS_NO_SW, RW> data0Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    HField& data0 = data0Accessor.getHostRW();
    data0 = 0;
    ModelArrayAccessor<CouplingFields::Q_SS_SW, RW> data1Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    data1Accessor.getHostRW() = data0;
    ModelArrayAccessor<CouplingFields::FWFLUX, RW> data2Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    data2Accessor.getHostRW() = data0;
    ModelArrayAccessor<CouplingFields::SFLUX, RW> data3Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    data3Accessor.getHostRW() = data0;*/

    // External SS* data
    ModelArrayAccessor<Protected::EXT_SSS, RW> sssExtAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sssExt = sssExtAccessor.getHostRW();
    sssExt = sss;

    ModelArrayAccessor<Protected::EXT_SST, RW> sstExtAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sstExt = sstExtAccessor.getHostRW();
    sstExt = sst + tOffset;

    SlabOcean slabOcean(couplingArrays);
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayAccessor<Protected::SLAB_QDW> qdwAccessor(ModelComponent::getStore());
    const HField& qdw = qdwAccessor.getHostRO();

    FloatType prec = 1e-8;
    REQUIRE(qdw[0]
        == doctest::Approx(tOffset * cpml[0] / SlabOcean::defaultRelaxationTime).epsilon(prec));

    ModelArrayAccessor<Protected::SLAB_SST> sstSlabAccessor(ModelComponent::getStore());
    // scope needed because we have to access sstSlab again after update
    {
        const HField& sstSlab = sstSlabAccessor.getHostRO();

        REQUIRE(sstSlab[0] != doctest::Approx(sst[0]).epsilon(prec / dt));
        REQUIRE(sstSlab[0] == doctest::Approx(sst[0] + dt * qdw[0] / cpml[0]).epsilon(prec));
    }

    ModelArrayAccessor<CouplingFields::Q_SS_SW, RW> qswNetAccessor(
        couplingArrays, RW, ModelArray::Type::H);
    HField& qswNet = qswNetAccessor.getHostRW();
    qswNet[0] = 15;
    ModelArrayAccessor<CouplingFields::Q_SS_NO_SW, RW> qNoSunAccessor(
        couplingArrays, RW, ModelArray::Type::H);
    HField& qNoSun = qNoSunAccessor.getHostRW();
    qNoSun[0] = -17.5;

    // Should not need to update anything else, as the slabOcean update only changes SLAB_SST
    slabOcean.update(tst);
    const HField& sstSlab = sstSlabAccessor.getHostRO();
    REQUIRE(sstSlab[0]
        == doctest::Approx(sst[0] - dt * (qswNet[0] + qNoSun[0] - qdw[0]) / cpml[0]).epsilon(prec));
}

TEST_CASE("Test Fdw")
{
    // 1000 s timestep
    TimestepTime tst = { TimePoint(), Duration("P0-0T0:16:40") };
    FloatType dt = tst.step.seconds();
    REQUIRE(dt == 1000);

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    ModelArrayStore couplingArrays;

    FloatType sOffset = 0.1;
    // Supply the data to the slab ocean
    ModelArrayAccessor<Protected::SSS, RW> sssAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sss = sssAccessor.getHostRW();
    sss = 32.;

    ModelArrayAccessor<Protected::SST, RW> sstAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sst = sstAccessor.getHostRW();
    sst = LinearFreezing()(sss[0]);

    ModelArrayAccessor<Protected::MLD, RW> mldAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& mld = mldAccessor.getHostRW();
    mld = 6.48;

    ModelArrayAccessor<Protected::ML_BULK_CP, RW> cpmlAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& cpml = cpmlAccessor.getHostRW();
    cpml = Water::cp * Water::rho * mld;

    /*
    ModelArrayAccessor<CouplingFields::Q_SS_NO_SW, RW> data0Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    HField& data0 = data0Accessor.getHostRW();
    data0 = 0;
    ModelArrayAccessor<CouplingFields::Q_SS_SW, RW> data1Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    data1Accessor.getHostRW() = data0;
    ModelArrayAccessor<CouplingFields::FWFLUX, RW> data2Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    data2Accessor.getHostRW() = data0;
    ModelArrayAccessor<CouplingFields::SFLUX, RW> data3Accessor(
        couplingArrays, RO, ModelArray::Type::H);
    data3Accessor.getHostRW() = data0;*/

    // External SS* data
    ModelArrayAccessor<Protected::EXT_SSS, RW> sssExtAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sssExt = sssExtAccessor.getHostRW();
    sssExt = sss + sOffset;

    ModelArrayAccessor<Protected::EXT_SST, RW> sstExtAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& sstExt = sstExtAccessor.getHostRW();
    sstExt = sst;

    SlabOcean slabOcean(couplingArrays);
    slabOcean.configure();
    slabOcean.update(tst);

    ModelArrayAccessor<Protected::SLAB_FDW> fdwAccessor(ModelComponent::getStore());
    const HField& fdw = fdwAccessor.getHostRO();

    FloatType prec = 1e-6;
    REQUIRE(fdw[0]
        == doctest::Approx(
            -sOffset / sss[0] * mld[0] * Water::rho / SlabOcean::defaultRelaxationTime)
               .epsilon(prec));
    // Test that the finiteelement.cpp calculation of fdw is not being used
    FloatType delS = -sOffset;
    FloatType timeS = SlabOcean::defaultRelaxationTime;
    FloatType ddt = tst.step.seconds();
    FloatType oldFdw = delS * mld[0] * Water::rho / (timeS * sss[0] - ddt * delS);
    REQUIRE(fdw[0] != doctest::Approx(oldFdw).epsilon(prec * 1e-6));

    ModelArrayAccessor<Protected::SLAB_SSS> sssSlabAccessor(ModelComponent::getStore());
    const HField& sssSlab = sssSlabAccessor.getHostRO();

    REQUIRE(sssSlab[0] != doctest::Approx(sss[0]).epsilon(prec / dt));
    REQUIRE(sssSlab[0]
        == doctest::Approx(sss[0] - (fdw[0] * dt) / (mld[0] * Water::rho + fdw[0] * dt))
               .epsilon(prec));

    ModelArrayAccessor<CouplingFields::FWFLUX, RW> snowMeltFluxAccessor(
        couplingArrays, RW, ModelArray::Type::H);
    HField& snowMeltFlux = snowMeltFluxAccessor.getHostRW();
    FloatType snowMelt = -1e-4;
    FloatType snowMeltVol = snowMelt * Ice::rhoSnow;
    snowMeltFlux = snowMeltVol / dt;
    slabOcean.update(tst);
    REQUIRE(sssSlab[0]
        == doctest::Approx(sss[0]
            + (snowMeltVol - fdw[0] * dt) / (mld[0] * Water::rho - snowMeltVol + fdw[0] * dt))
               .epsilon(prec));
}
TEST_SUITE_END();
} /* namespace Nextsim */
