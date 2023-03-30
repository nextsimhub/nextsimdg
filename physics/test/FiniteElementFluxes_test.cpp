/*!
 * @file FiniteElementFluxes_test.cpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/FiniteElementFluxes.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Module.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_CASE("Melting conditions", "[FiniteElementFluxes]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;
    config << "Nextsim::IIceAlbedo = Nextsim::CCSMIceAlbedo" << std::endl;
    config << std::endl;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1., 32., 10.25);
    ocnBdy.setData(ModelState().data);

    class AtmosphereData : public ModelComponent {
    public:
        AtmosphereData() {
            getStore().registerArray(Protected::T_AIR, &tair);
            getStore().registerArray(Protected::DEW_2M, &tdew);
            getStore().registerArray(Protected::P_AIR, &pair);
            getStore().registerArray(Protected::WIND_SPEED, &windSpeed);
            getStore().registerArray(Protected::SW_IN, &sw_in);
            getStore().registerArray(Protected::LW_IN, &lw_in);
        }
        void setData(const ModelState::DataMap& state) override
        {
            tair.resize();
            tdew.resize();
            pair.resize();
            windSpeed.resize();
            sw_in.resize();
            lw_in.resize();

            tair = 3;
            tdew = 2;
            pair = 100000.;
            windSpeed = 5;
            sw_in = 50;
            lw_in = 330;
        }
        std::string getName() const override { return "AtmData"; }
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

    private:
        HField tair;
        HField tdew;
        HField pair;
        HField windSpeed;
        HField sw_in;
        HField lw_in;
        HField snowfall;
    } atmState;
    atmState.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            getStore().registerArray(Protected::H_ICE, &hice);
            getStore().registerArray(Protected::C_ICE, &cice);
            getStore().registerArray(Protected::H_SNOW, &hsnow);
            getStore().registerArray(Protected::T_ICE, &tice0);
            getStore().registerArray(Protected::HTRUE_ICE, &hice0);
            getStore().registerArray(Protected::HTRUE_SNOW, &hsnow0);

        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -1.;

            hice0[0] = hice[0] / cice[0];
            hsnow0[0] = hsnow[0] / cice[0];
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        HField hice0;  // ice averaged ice thickness
        HField hsnow0; // ice averaged snow thickness
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    HField qow;
    qow.resize();
    ModelComponent::getStore().registerArray(Shared::Q_OW, &qow, RW);

    HField qia;
    qia.resize();
    ModelComponent::getStore().registerArray(Shared::Q_IA, &qia, RW);

    HField dqia_dt;
    dqia_dt.resize();
    ModelComponent::getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);

    HField subl;
    subl.resize();
    ModelComponent::getStore().registerArray(Shared::SUBLIM, &subl, RW);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    // OceanState is independently updated
    FiniteElementFluxes fef;
    fef.configure();
    fef.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    fef.update(tst);

    double prec = 1e-5;
    REQUIRE(qow[0] == Approx(-109.923).epsilon(prec));
    REQUIRE(qia[0] == Approx(-84.5952).epsilon(prec));
    REQUIRE(dqia_dt[0] == Approx(19.7016).epsilon(prec));
    REQUIRE(subl[0] == Approx(-7.3858e-06).epsilon(prec));
}

TEST_CASE("Freezing conditions", "[FiniteElementFluxes]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;
    config << "Nextsim::IIceAlbedo = Nextsim::CCSMIceAlbedo" << std::endl;
    config << std::endl;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1.75, 32., 10.25);
    ocnBdy.setData(ModelState().data);

    class AtmosphereData : public ModelComponent {
    public:
        AtmosphereData()
        {
            getStore().registerArray(Protected::T_AIR, &tair);
            getStore().registerArray(Protected::DEW_2M, &tdew);
            getStore().registerArray(Protected::P_AIR, &pair);
            getStore().registerArray(Protected::WIND_SPEED, &windSpeed);
            getStore().registerArray(Protected::SW_IN, &sw_in);
            getStore().registerArray(Protected::LW_IN, &lw_in);
        }
        void setData(const ModelState::DataMap& state) override
        {
            tair.resize();
            tdew.resize();
            pair.resize();
            windSpeed.resize();
            sw_in.resize();
            lw_in.resize();
            tair = -12;
            tdew = -12;
            pair = 100000.;
            windSpeed = 5;
            sw_in = 0;
            lw_in = 265;
        }
        std::string getName() const override { return "AtmData"; }
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

    private:
        HField tair;
        HField tdew;
        HField pair;
        HField windSpeed;
        HField sw_in;
        HField lw_in;
        HField snowfall;
    } atmState;
    atmState.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            getStore().registerArray(Protected::H_ICE, &hice);
            getStore().registerArray(Protected::C_ICE, &cice);
            getStore().registerArray(Protected::H_SNOW, &hsnow);
            getStore().registerArray(Protected::T_ICE, &tice0);
            getStore().registerArray(Protected::HTRUE_ICE, &hice0);
            getStore().registerArray(Protected::HTRUE_SNOW, &hsnow0);
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -9.;

            hice0[0] = hice[0] / cice[0];
            hsnow0[0] = hsnow[0] / cice[0];
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        HField hice0;  // ice averaged ice thickness
        HField hsnow0; // ice averaged snow thickness

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    HField qow;
    qow.resize();
    ModelComponent::getStore().registerArray(Shared::Q_OW, &qow, RW);

    HField qia;
    qia.resize();
    ModelComponent::getStore().registerArray(Shared::Q_IA, &qia, RW);

    HField dqia_dt;
    dqia_dt.resize();
    ModelComponent::getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);

    HField subl;
    subl.resize();
    ModelComponent::getStore().registerArray(Shared::SUBLIM, &subl, RW);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    // OceanState is independently updated
    FiniteElementFluxes fef;
    fef.configure();
    ocnBdy.updateBefore(tst);
    fef.setData(ModelState().data);
    fef.update(tst);

    double prec = 1e-5;
    REQUIRE(qow[0] == Approx(143.266).epsilon(prec));
    REQUIRE(qia[0] == Approx(42.2955).epsilon(prec));
    REQUIRE(dqia_dt[0] == Approx(16.7615).epsilon(prec));
    REQUIRE(subl[0] == Approx(2.15132e-6).epsilon(prec));
}
}
