/*!
 * @file ThermoWintonTemperature_test.cpp
 *
 * @date Nov 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/IFluxCalculation.hpp"
#include "include/ThermoWinton.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/constants.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/IAtmosphereBoundary.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"
#include "include/UniformOcean.hpp"

namespace Nextsim {

TEST_CASE("Melting conditions", "[ThermoWinton]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 3 });

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

    ThermoWinton twin;
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : tice0(ModelArray::Type::Z)
            , tice(getStore())
        {
            getStore().registerArray(Protected::HTRUE_ICE, &hice0);
            getStore().registerArray(Protected::C_ICE, &cice0);
            getStore().registerArray(Protected::HTRUE_SNOW, &hsnow0);
            getStore().registerArray(Protected::SW_IN, &sw_in);
            getStore().registerArray(Protected::T_ICE, &tice0);

            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice0[0] = 0.5;
            hice0[0] = 0.1 / cice0[0]; // Here we are using the true thicknesses
            hsnow0[0] = 0.01 / cice0[0];
            sw_in[0] = -10.1675; // Net shortwave flux from incident 50 W/m^2
            tice0[0] = -1;
            tice0[1] = -1;
            tice0[2] = -1;
            tice.data().setData(tice0);

            hice = hice0;
            cice = cice0;
            hsnow = hsnow0;

        }

        HField hice0;
        HField cice0;
        HField hsnow0;
        HField sw_in;
        ZField tice0;
        ModelArrayRef<Shared::T_ICE, RW> tice; // From IIceThermodynamics

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } initCond;
    initCond.setData(ModelState().data);

    UniformOcean oceanData(-1, 32., 4.29151e7/(Water::rho * Water::cp));
    oceanData.setQio(53717.8);
    oceanData.setData(ModelState().data);

    class AtmosphereState : public IAtmosphereBoundary {
    public:
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            snow[0] = 0.00;
            qow[0] = -109.923;
            qia[0] = -84.5952;
            dqia_dt[0] = 19.7016;
            subl[0] = -7.3858e-06;
        }
    } atmosState;
    atmosState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    twin.configure();
    twin.update(tst);
    ModelArrayRef<Shared::T_ICE, RO> tice(ModelComponent::getStore());
    ModelArrayRef<Shared::Q_IC, RO> qic(ModelComponent::getStore());

    double prec = 1e-5;

    REQUIRE(tice[0] == Approx(0.0).margin(prec));
    REQUIRE(tice[1] == Approx(-0.999453).margin(prec));
    REQUIRE(tice[2] == Approx(-0.275).epsilon(prec));
    //    REQUIRE(qic[0] == Approx(-4.60879).epsilon(prec));
}

TEST_CASE("Freezing conditions", "[ThermoWinton]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 3 });

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

    ThermoWinton twin;
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : tice0(ModelArray::Type::Z)
        , tice(getStore())
        {
            getStore().registerArray(Protected::HTRUE_ICE, &hice0);
            getStore().registerArray(Protected::C_ICE, &cice0);
            getStore().registerArray(Protected::HTRUE_SNOW, &hsnow0);
            getStore().registerArray(Protected::SNOW, &snow);
            getStore().registerArray(Protected::SW_IN, &sw_in);
            getStore().registerArray(Protected::T_ICE, &tice0);

            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);

        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice0[0] = 0.5;
            hice0[0] = 0.1 / cice0[0]; // Here we are using the true thicknesses
            hsnow0[0] = 0.01 / cice0[0];
            snow[0] = 1e-3;
            sw_in[0] = 0;
            tice0[0] = -9.;
            tice0[1] = -9.;
            tice0[2] = -9.;
            tice.data().setData(tice0);

            hice = hice0;
            cice = cice0;
            hsnow = hsnow0;

        }

        HField hice0;
        HField cice0;
        HField hsnow0;
        HField snow;
        HField sw_in;
        ZField tice0;
        ModelArrayRef<Shared::T_ICE, RW> tice; // From IIceThermodynamics

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState().data);

    UniformOcean oceanData(-1.75, 32., 4.29151e7/(Water::rho * Water::cp));
    oceanData.setQio(73.9465);
    oceanData.setData(ModelState().data);

    class AtmosphereState : public IAtmosphereBoundary {
    public:
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            snow[0] = 0.00;
            qow[0] = 143.266;
            qia[0] = 42.2955;
            dqia_dt[0] = 16.7615;
            subl[0] = 2.15132e-6;
        }
    } atmosState;
    atmosState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    twin.configure();
    twin.update(tst);

    ModelArrayRef<Shared::T_ICE, RO> tice(ModelComponent::getStore());
    ModelArrayRef<Shared::Q_IC, RO> qic(ModelComponent::getStore());

    double prec = 1e-5;

    REQUIRE(tice[0] == Approx(-10.5129).epsilon(prec));
    REQUIRE(tice[1] == Approx(-9.00726).epsilon(prec));
    REQUIRE(tice[2] == Approx(-8.20454).epsilon(prec));
//    REQUIRE(qic[0] == Approx(44.4839).epsilon(prec));
}
}
