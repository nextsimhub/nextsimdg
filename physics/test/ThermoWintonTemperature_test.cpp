/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
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

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ThermoWinton twin;

#define TICE -1. // Cannot reference local var in the function definition

    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : sw_inAccessor(getStore(), RO)
            , hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            ciceAccessor.getHostRW()[0] = 0.5;
            hiceAccessor.getHostRW()[0] = 0.1;
            hsnowAccessor.getHostRW()[0] = 0.01;
            sw_inAccessor.getHostRW()[0] = -10.1675; // Net shortwave flux from incident 50 W/m^2
        }

        ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
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
            snowAccessor.getHostRW()[0] = 0.00;
            qowAccessor.getHostRW()[0] = -109.923;
            qiaAccessor.getHostRW()[0] = -85.6364;
            dqia_dtAccessor.getHostRW()[0] = 19.7016;
            sublAccessor.getHostRW()[0] = -7.3858e-06;
            penSWAccessor.getHostRW()[0] = 1.04125;
        }
    } atmosState;
    atmosState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    // Data for the Winton ice temperatures
    ModelState::DataMap iceTempState = {
        { tsurfName, ModelArray(ModelArray::Type::H) },
        { ThermoWinton::tInteriorName, ModelArray(ModelArray::Type::H) },
        { ThermoWinton::tBottomName, ModelArray(ModelArray::Type::H) },
    };
    iceTempState.at(tsurfName) = TICE;
    iceTempState.at(ThermoWinton::tInteriorName) = TICE;
    iceTempState.at(ThermoWinton::tBottomName) = TICE;
#undef TICE

    twin.configure();
    twin.setData(iceTempState);
    twin.update(tst);

    ModelState::DataMap output = twin.getStatePrognostic().data;
    REQUIRE(output.count(tsurfName) != 0);
    REQUIRE(output.count(ThermoWinton::tInteriorName) != 0);
    REQUIRE(output.count(ThermoWinton::tBottomName) != 0);

    double prec = 1e-5;

    REQUIRE(output.at(tsurfName)[0] == doctest::Approx(0.0).epsilon(prec));
    REQUIRE(output.at(ThermoWinton::tInteriorName)[0] == doctest::Approx(-0.999261).epsilon(prec));
    REQUIRE(output.at(ThermoWinton::tBottomName)[0] == doctest::Approx(-0.275).epsilon(prec));
}

TEST_CASE("Freezing conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ThermoWinton twin;
#define TICE -9.
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : snowAccessor(getStore(), RO)
            , sw_inAccessor(getStore(), RO)
            , hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice0[0] = 0.5;
            hice0[0] = 0.1;
            hsnow0[0] = 0.01;
            snowAccessor.getHostRW()[0] = 1e-3;
            sw_inAccessor.getHostRW()[0] = 0;

            hiceAccessor.getHostRW() = hice0;
            ciceAccessor.getHostRW() = cice0;
            hsnowAccessor.getHostRW() = hsnow0;
        }

        HField hice0;
        HField cice0;
        HField hsnow0;
        ModelArrayAccessor<Protected::SNOW, RW> snowAccessor;
        ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
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
            snowAccessor.getHostRW()[0] = 0.00;
            qowAccessor.getHostRW()[0] = 143.266;
            qiaAccessor.getHostRW()[0] = 42.2955;
            dqia_dtAccessor.getHostRW()[0] = 16.7615;
            sublAccessor.getHostRW()[0] = 2.15132e-6;
            penSWAccessor.getHostRW()[0] = 0.;
        }
    } atmosState;
    atmosState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    // Data for the Winton ice temperatures
    ModelState::DataMap iceTempState = {
        { tsurfName, ModelArray(ModelArray::Type::H) },
        { ThermoWinton::tInteriorName, ModelArray(ModelArray::Type::H) },
        { ThermoWinton::tBottomName, ModelArray(ModelArray::Type::H) },
    };
    iceTempState.at(tsurfName) = TICE;
    iceTempState.at(ThermoWinton::tInteriorName) = TICE;
    iceTempState.at(ThermoWinton::tBottomName) = TICE;
#undef TICE

    twin.configure();
    twin.setData(iceTempState);
    twin.update(tst);

    ModelState::DataMap output = twin.getStatePrognostic().data;
    REQUIRE(output.count(tsurfName) != 0);
    REQUIRE(output.count(ThermoWinton::tInteriorName) != 0);
    REQUIRE(output.count(ThermoWinton::tBottomName) != 0);

    double prec = 1e-5;

    REQUIRE(output.at(tsurfName)[0] == doctest::Approx(-10.5129).epsilon(prec));
    REQUIRE(output.at(ThermoWinton::tInteriorName)[0] == doctest::Approx(-9.00726).epsilon(prec));
    REQUIRE(output.at(ThermoWinton::tBottomName)[0] == doctest::Approx(-8.20454).epsilon(prec));
    //    REQUIRE(qic[0] == doctest::Approx(44.4839).epsilon(prec));
}

TEST_CASE("No ice do nothing")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ThermoWinton twin;
#define TICE 0.
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : snowAccessor(getStore(), RO)
            , sw_inAccessor(getStore(), RO)
            , hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice0[0] = 0;
            hice0[0] = 0;
            hsnow0[0] = 0;
            snowAccessor.getHostRW()[0] = 0;
            sw_inAccessor.getHostRW()[0] = 0;

            hiceAccessor.getHostRW() = hice0;
            ciceAccessor.getHostRW() = cice0;
            hsnowAccessor.getHostRW() = hsnow0;
        }

        HField hice0;
        HField cice0;
        HField hsnow0;
        ModelArrayAccessor<Protected::SNOW, RW> snowAccessor;
        ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } atmoState;
    atmoState.setData(ModelState().data);

    class OceanState : public IOceanBoundary {
    public:
        void setData(const ModelState::DataMap& ms) override
        {
            IOceanBoundary::setData(ms);
            sstAccessor.getHostRW()[0] = 1.75;
            HField& sss = sssAccessor.getHostRW();
            sss[0] = 32.;
            tfAccessor.getHostRW()[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            cpmlAccessor.getHostRW()[0] = 4.29151e7;
            qioAccessor.getHostRW()[0] = 0;
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
            snowAccessor.getHostRW()[0] = 0.00;
            qowAccessor.getHostRW()[0] = 143.266;
            qiaAccessor.getHostRW()[0] = 42.2955;
            dqia_dtAccessor.getHostRW()[0] = 16.7615;
            sublAccessor.getHostRW()[0] = 2.15132e-6;
        }
    } atmosState;
    atmosState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    // Data for the Winton ice temperatures
    ModelState::DataMap iceTempState = {
        { tsurfName, ModelArray(ModelArray::Type::H) },
        { ThermoWinton::tInteriorName, ModelArray(ModelArray::Type::H) },
        { ThermoWinton::tBottomName, ModelArray(ModelArray::Type::H) },
    };
    iceTempState.at(tsurfName) = TICE;
    iceTempState.at(ThermoWinton::tInteriorName) = TICE;
    iceTempState.at(ThermoWinton::tBottomName) = TICE;
#undef TICE

    twin.configure();
    twin.update(tst);

    ModelArrayAccessor<Shared::H_ICE_DG, RO> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();

    //    double prec = 1e-5;

    REQUIRE(hice[0] == 0);
    REQUIRE(cice[0] == 0);
}

TEST_SUITE_END();

}
