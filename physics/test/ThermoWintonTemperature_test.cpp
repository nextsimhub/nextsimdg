/*!
 * @file ThermoWintonTemperature_test.cpp
 *
 * @date 29 Jan 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <sstream>

#include "include/IFluxCalculation.hpp"
#include "include/ThermoWinton.hpp"

#include "include/Configurator.hpp"
#include "include/IAtmosphereBoundary.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/IOceanBoundary.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ThermoWintonTemperature");
TEST_CASE("Melting conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 3 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ThermoWinton twin;
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : tice0(ModelArray::Type::Z)
        {
            getStore().registerArray(Protected::SW_IN, &sw_in, RO);
            getStore().registerArray(Shared::T_ICE, &tice0, RW);

            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice[0] = 0.5;
            hice[0] = 0.1;
            hsnow[0] = 0.01;
            sw_in[0] = -10.1675; // Net shortwave flux from incident 50 W/m^2
            tice0[0] = -1;
            tice0[1] = -1;
            tice0[2] = -1;
        }

        HField sw_in;
        ZField tice0;

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } initCond;
    initCond.setData(ModelState().data);

    UniformOcean oceanData(-1, 32., 4.29151e7 / (Water::rho * Water::cp));
    oceanData.setQio(53717.8);
    oceanData.setData(ModelState().data);

    class AtmosphereState : public IAtmosphereBoundary {
    public:
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            snow[0] = 0.00;
            qow[0] = -109.923;
            qia[0] = -85.6364;
            dqia_dt[0] = 19.7016;
            subl[0] = -7.3858e-06;
            penSW[0] = 1.04125;
        }
    } atmosState;
    atmosState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    twin.configure();
    twin.update(tst);
    ModelArrayRef<Shared::T_ICE, RO> tice(ModelComponent::getStore());
    ModelArrayRef<Shared::Q_IC, RO> qic(ModelComponent::getStore());

    double prec = 1e-5;

    REQUIRE(tice[0] == doctest::Approx(0.0).epsilon(prec));
    REQUIRE(tice[1] == doctest::Approx(-0.999261).epsilon(prec));
    REQUIRE(tice[2] == doctest::Approx(-0.275).epsilon(prec));
    //    REQUIRE(qic[0] == doctest::Approx(-4.60879).epsilon(prec));
}

TEST_CASE("Freezing conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 3 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ThermoWinton twin;
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : tice0(ModelArray::Type::Z)
        {
            getStore().registerArray(Protected::SNOW, &snow, RO);
            getStore().registerArray(Protected::SW_IN, &sw_in, RO);
            getStore().registerArray(Shared::T_ICE, &tice0, RW);

            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice[0] = 0.5;
            hice[0] = 0.1;
            hsnow[0] = 0.01;
            snow[0] = 1e-3;
            sw_in[0] = 0;
            tice0[0] = -9.;
            tice0[1] = -9.;
            tice0[2] = -9.;
        }

        HField snow;
        HField sw_in;
        ZField tice0;

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState().data);

    UniformOcean oceanData(-1.75, 32., 4.29151e7 / (Water::rho * Water::cp));
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

    REQUIRE(tice[0] == doctest::Approx(-10.5129).epsilon(prec));
    REQUIRE(tice[1] == doctest::Approx(-9.00726).epsilon(prec));
    REQUIRE(tice[2] == doctest::Approx(-8.20454).epsilon(prec));
    //    REQUIRE(qic[0] == doctest::Approx(44.4839).epsilon(prec));
}

TEST_CASE("No ice do nothing")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 3 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ThermoWinton twin;
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : tice0(ModelArray::Type::Z)
        {
            getStore().registerArray(Protected::SNOW, &snow, RO);
            getStore().registerArray(Protected::SW_IN, &sw_in, RO);
            getStore().registerArray(Shared::T_ICE, &tice0, RW);

            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice[0] = 0;
            hice[0] = 0;
            hsnow[0] = 0;
            snow[0] = 0;
            sw_in[0] = 0;
            tice0[0] = 0;
            tice0[1] = 0;
            tice0[2] = 0.;
        }

        HField snow;
        HField sw_in;
        ZField tice0;

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState().data);

    class OceanState : public IOceanBoundary {
    public:
        void setData(const ModelState::DataMap& ms) override
        {
            IOceanBoundary::setData(ms);
            sst[0] = 1.75;
            sss[0] = 32.;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            cpml[0] = 4.29151e7;
            qio[0] = 0;
        }
        void updateBefore(const TimestepTime& tst) override { }
        void updateAfter(const TimestepTime& tst) override { }
    } oceanData;
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

    ModelArrayRef<Shared::H_ICE, RO> hice(ModelComponent::getStore());
    ModelArrayRef<Shared::C_ICE, RO> cice(ModelComponent::getStore());

    double prec = 1e-5;

    REQUIRE(hice[0] == 0);
    REQUIRE(cice[0] == 0);
}

TEST_SUITE_END();

}
