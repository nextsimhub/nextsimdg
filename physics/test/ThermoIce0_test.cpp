/*!
 * @file ThermoIce0Temperature_test.cpp
 *
 * @date 19 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <sstream>

#include "include/ThermoIce0.hpp"
#include "include/IFluxCalculation.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/FreezingPointModule.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ThermoIce0");
/*
 * Test that ice below the minimum ice threshold is eliminated.
 */
TEST_CASE("Threshold ice")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    // Class derived from ModelComponent providing the physical data for the test
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RO);
            getStore().registerArray(Protected::SST, &sst, RO);
            getStore().registerArray(Protected::SSS, &sss, RO);
            getStore().registerArray(Protected::TF, &tf, RO);
            getStore().registerArray(Protected::SNOW, &snow, RO);
            getStore().registerArray(Protected::ML_BULK_CP, &mlbhc, RO);
            getStore().registerArray(Shared::Q_IO, &qio, RW);
            getStore().registerArray(Shared::Q_OW, &qow, RW);
            getStore().registerArray(Shared::Q_IA, &qia, RW);
            getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);
            getStore().registerArray(Shared::SUBLIM, &subl, RW);
            getStore().registerArray(Shared::Q_PEN_SW, &penSW, RW);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice[0] = 0.99;
            hice[0] = 0.001 / cice[0]; // Here we are using the true thicknesses
            hsnow[0] = 0.;
            sss[0] = 32.;
            sst[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            snow[0] = 0.;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;
            qio[0] = 0.;
            qow[0] = 0;
            qia[0] = 0;
            dqia_dt[0] = 0;
            subl[0] = 0;
            penSW[0] = 0;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity
        HField qio;
        HField qow;
        HField qia;
        HField dqia_dt;
        HField subl;
        HField penSW;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState::DataMap());

    // Supply the atmospheric boundary arrays, without an entire
    // IAtmosphereBoundary implementation
    HField qow;
    qow.resize();
    ModelComponent::getStore().registerArray(Shared::Q_OW, &qow, RW);

    HField qia;
    qia.resize();
    ModelComponent::getStore().registerArray(Shared::Q_IA, &qia, RW);

    HField dqia_dt;
    dqia_dt.resize();
    ModelComponent::getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);

    HField qic;
    qic.resize();
    ModelComponent::getStore().registerArray(Shared::Q_IC, &qic, RW);

    HField subl;
    subl.resize();
    ModelComponent::getStore().registerArray(Shared::SUBLIM, &subl, RW);

    // An implementation of IFluxCalculation that returns zero fluxes
    class FluxData : public IFluxCalculation {
    public:
        FluxData()
            : IFluxCalculation()
        {
        }
        std::string getName() const override { return "FluxData"; }

        void setData(const ModelState::DataMap&) override
        {
            qow[0] = 0;
            qia[0] = 0;
            dqia_dt[0] = 0;
            subl[0] = 0;
            penSW[0] = 0;
        }

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

        void update(const TimestepTime&) override { }
    } fluxData;

    fluxData.setData(ModelState::DataMap());

    TimestepTime tst = { TimePoint("2000-01-01T00:00:00"), Duration(600) };
    ThermoIce0 ti0t;
    ti0t.configure();
    ti0t.setData(ModelState::DataMap());
    ti0t.update(tst);

    ModelArrayRef<Shared::H_ICE> hice(ModelComponent::getStore());
    // So little ice should be reduced to zero
    REQUIRE(hice[0] == 0.);
    ModelArrayRef<Shared::C_ICE> cice(ModelComponent::getStore());
    REQUIRE(cice[0] == 0.);

}
TEST_SUITE_END();
}
