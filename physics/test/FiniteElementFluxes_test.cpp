/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include <doctest/doctest.h>
#include <sstream>

#include "include/FiniteElementFluxes.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("FiniteElementFluxes");
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

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1., 32., 10.25);
    ocnBdy.setData(ModelState().data);

    class AtmosphereData : public ModelComponent {
    public:
        AtmosphereData()
            : tairAccessor(getStore(), RO)
            , tdewAccessor(getStore(), RO)
            , pairAccessor(getStore(), RO)
            , windSpeedAccessor(getStore(), RO)
            , u_airAccessor(getStore(), RO)
            , v_airAccessor(getStore(), RO)
            , sw_inAccessor(getStore(), RO)
            , lw_inAccessor(getStore(), RO)
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            HField& tair = tairAccessor.getHostRW();
            tair.resize();
            HField& tdew = tdewAccessor.getHostRW();
            tdew.resize();
            HField& pair = pairAccessor.getHostRW();
            pair.resize();
            HField& windSpeed = windSpeedAccessor.getHostRW();
            windSpeed.resize();
            HField& u_air = u_airAccessor.getHostRW();
            u_air.resize();
            HField& v_air = v_airAccessor.getHostRW();
            v_air.resize();
            HField& sw_in = sw_inAccessor.getHostRW();
            sw_in.resize();
            HField& lw_in = lw_inAccessor.getHostRW();
            lw_in.resize();

            tair = 3;
            tdew = 2;
            pair = 100000.;
            windSpeed = 5;
            u_air = 3;
            v_air = 4;
            sw_in = 50;
            lw_in = 330;
        }
        std::string getName() const override { return "AtmData"; }

    private:
        ModelArrayAccessor<Protected::T_AIR, RW> tairAccessor;
        ModelArrayAccessor<Protected::DEW_2M, RW> tdewAccessor;
        ModelArrayAccessor<Protected::P_AIR, RW> pairAccessor;
        ModelArrayAccessor<Protected::WIND_SPEED, RW> windSpeedAccessor;
        ModelArrayAccessor<Protected::WIND_U, RW> u_airAccessor;
        ModelArrayAccessor<Protected::WIND_V, RW> v_airAccessor;
        ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;
        ModelArrayAccessor<Protected::LW_IN, RW> lw_inAccessor;
        HField snowfall;
    } atmState;
    atmState.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
            : hiceAccessor(getStore(), RO)
            , ciceAccessor(getStore(), RO)
            , hsnowAccessor(getStore(), RO)
            , tsurfAccessor(getStore(), RO)
        {
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW()[0] = 0.5;
            hiceAccessor.getHostRW()[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnowAccessor.getHostRW()[0] = 0.01;
            tsurfAccessor.getHostRW()[0] = -1.;
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
        ModelArrayAccessor<Protected::T_SURF, RW> tsurfAccessor;
    } iceState;
    iceState.setData(ModelState().data);

    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor(ModelComponent::getStore(), RW);
    qowAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::Q_IA, RW> qiaAccessor(ModelComponent::getStore(), RW);
    qiaAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::Q_PEN_SW, RW> penSWAccessor(ModelComponent::getStore(), RW);
    penSWAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::Q_SW_OW, RW> qsw_owAccessor(ModelComponent::getStore(), RW);
    qsw_owAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::Q_SW_BASE, RW> qsw_baseAccessor(ModelComponent::getStore(), RW);
    qsw_baseAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::DQIA_DT, RW> dqia_dtAccessor(ModelComponent::getStore(), RW);
    dqia_dtAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::SUBLIM, RW> sublAccessor(ModelComponent::getStore(), RW);
    sublAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::EVAP, RW> evapAccessor(ModelComponent::getStore(), RW);
    evapAccessor.getHostRW().resize();

    ModelArrayAccessor<Shared::OW_STRESS_X, RW> tauXAccessor(ModelComponent::getStore(), RW);
    tauXAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::OW_STRESS_Y, RW> tauYAccessor(ModelComponent::getStore(), RW);
    tauYAccessor.getHostRW().resize();


    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    // OceanState is independently updated
    FiniteElementFluxes fef;
    fef.configure();
    fef.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    fef.update(tst);

    double prec = 1e-5;
    REQUIRE(qowAccessor.getHostRO()[0] == doctest::Approx(-109.923).epsilon(prec));
    REQUIRE(qiaAccessor.getHostRO()[0] == doctest::Approx(-85.6364).epsilon(prec));
    REQUIRE(dqia_dtAccessor.getHostRO()[0] == doctest::Approx(19.7016).epsilon(prec));
    REQUIRE(sublAccessor.getHostRO()[0] == doctest::Approx(-7.3858e-06).epsilon(prec));
    REQUIRE(tauXAccessor.getHostRO()[0] == doctest::Approx(1.89732e-2).epsilon(prec));
    REQUIRE(tauYAccessor.getHostRO()[0] == doctest::Approx(2.52976e-2).epsilon(prec));
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

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1.75, 32., 10.25);
    ocnBdy.setData(ModelState().data);

    class AtmosphereData : public ModelComponent {
    public:
        AtmosphereData()
            : tairAccessor(getStore(), RO)
            , tdewAccessor(getStore(), RO)
            , pairAccessor(getStore(), RO)
            , windSpeedAccessor(getStore(), RO)
            , u_airAccessor(getStore(), RO)
            , v_airAccessor(getStore(), RO)
            , sw_inAccessor(getStore(), RO)
            , lw_inAccessor(getStore(), RO)
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            HField& tair = tairAccessor.getHostRW();
            tair.resize();
            HField& tdew = tdewAccessor.getHostRW();
            tdew.resize();
            HField& pair = pairAccessor.getHostRW();
            pair.resize();
            HField& windSpeed = windSpeedAccessor.getHostRW();
            windSpeed.resize();
            HField& u_air = u_airAccessor.getHostRW();
            u_air.resize();
            HField& v_air = v_airAccessor.getHostRW();
            v_air.resize();
            HField& sw_in = sw_inAccessor.getHostRW();
            sw_in.resize();
            HField& lw_in = lw_inAccessor.getHostRW();
            lw_in.resize();
            tair = -12;
            tdew = -12;
            pair = 100000.;
            windSpeed = 5;
            u_air = 3;
            v_air = 4;
            sw_in = 0;
            lw_in = 265;
        }
        std::string getName() const override { return "AtmData"; }

    private:
        ModelArrayAccessor<Protected::T_AIR, RW> tairAccessor;
        ModelArrayAccessor<Protected::DEW_2M, RW> tdewAccessor;
        ModelArrayAccessor<Protected::P_AIR, RW> pairAccessor;
        ModelArrayAccessor<Protected::WIND_SPEED, RW> windSpeedAccessor;
        ModelArrayAccessor<Protected::WIND_U, RW> u_airAccessor;
        ModelArrayAccessor<Protected::WIND_V, RW> v_airAccessor;
        ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;
        ModelArrayAccessor<Protected::LW_IN, RW> lw_inAccessor;
        HField snowfall;
    } atmState;
    atmState.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
            : hiceAccessor(getStore(), RO)
            , ciceAccessor(getStore(), RO)
            , hsnowAccessor(getStore(), RO)
            , tsurfAccessor(getStore(), RO)
        {
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW()[0] = 0.5;
            hiceAccessor.getHostRW()[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnowAccessor.getHostRW()[0] = 0.01;
            tsurfAccessor.getHostRW()[0] = -9.;
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
        ModelArrayAccessor<Protected::T_SURF, RW> tsurfAccessor;

    } iceState;
    iceState.setData(ModelState().data);

    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor(ModelComponent::getStore(), RW);
    qowAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::Q_IA, RW> qiaAccessor(ModelComponent::getStore(), RW);
    qiaAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::Q_PEN_SW, RW> penSWAccessor(ModelComponent::getStore(), RW);
    penSWAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::DQIA_DT, RW> dqia_dtAccessor(ModelComponent::getStore(), RW);
    dqia_dtAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::SUBLIM, RW> sublAccessor(ModelComponent::getStore(), RW);
    sublAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::EVAP, RW> evapAccessor(ModelComponent::getStore(), RW);
    evapAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::OW_STRESS_X, RW> tauXAccessor(ModelComponent::getStore(), RW);
    tauXAccessor.getHostRW().resize();
    ModelArrayAccessor<Shared::OW_STRESS_Y, RW> tauYAccessor(ModelComponent::getStore(), RW);
    tauYAccessor.getHostRW().resize();

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    // OceanState is independently updated
    FiniteElementFluxes fef;
    fef.configure();
    ocnBdy.updateBefore(tst);
    fef.setData(ModelState().data);
    fef.update(tst);

    double prec = 1e-5;
    REQUIRE(qowAccessor.getHostRO()[0] == doctest::Approx(143.266).epsilon(prec));
    REQUIRE(qiaAccessor.getHostRO()[0] == doctest::Approx(42.2955).epsilon(prec));
    REQUIRE(dqia_dtAccessor.getHostRO()[0] == doctest::Approx(16.7615).epsilon(prec));
    REQUIRE(sublAccessor.getHostRO()[0] == doctest::Approx(2.15132e-6).epsilon(prec));
    REQUIRE(tauXAccessor.getHostRO()[0] == doctest::Approx(2.00279e-2).epsilon(prec));
    REQUIRE(tauYAccessor.getHostRO()[0] == doctest::Approx(2.67038e-2).epsilon(prec));
}
TEST_SUITE_END();

}
