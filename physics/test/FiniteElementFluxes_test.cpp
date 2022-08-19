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
#include "include/IFreezingPointModule.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

#include "include/OceanStateModule.hpp"

namespace Nextsim {

TEST_CASE("Melting conditions", "[FiniteElementFluxes]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

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

    class OceanData : public OceanState {
    public:
        OceanData()
            : OceanState()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            sst[0] = -1.;
            sss[0] = 32.;
            mld[0] = 10.25;
        }
        void updateSpecial(const TimestepTime& tst) override { }
    };
    Module::Module<OceanState>::setExternalImplementation(Module::newImpl<OceanState, OceanData>);

    class AtmosphereData : public AtmosphereState {
    public:
        AtmosphereData()
            : AtmosphereState()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            tair[0] = 3;
            tdew[0] = 2;
            pair[0] = 100000.;
            windSpeed[0] = 5;
            sw_in[0] = 50;
            lw_in[0] = 330;
            snowfall[0] = 0;
        }
        void updateSpecial(const TimestepTime& tst) override { }
    };
    Module::Module<AtmosphereState>::setExternalImplementation(
        Module::newImpl<AtmosphereState, AtmosphereData>);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -1.;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    FiniteElementFluxCalc fefc;
    fefc.configure();
    fefc.setData(ModelState().data);
    fefc.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::Q_OW, RO> qow;
    ModelArrayRef<ModelComponent::SharedArray::Q_IO, RO> qio;
    ModelArrayRef<ModelComponent::SharedArray::Q_IA, RO> qia;
    ModelArrayRef<ModelComponent::SharedArray::DQIA_DT, RO> dqia_dt;
    ModelArrayRef<ModelComponent::SharedArray::SUBLIM, RO> subl;

    double prec = 1e-5;
    REQUIRE(qow[0] == Approx(-109.923).epsilon(prec));
    REQUIRE(qio[0] == Approx(53717.8).epsilon(prec));
    REQUIRE(qia[0] == Approx(-84.5952).epsilon(prec));
    REQUIRE(dqia_dt[0] == Approx(19.7016).epsilon(prec));
    REQUIRE(subl[0] == Approx(-7.3858e-06).epsilon(prec));
}

TEST_CASE("Freezing conditions", "[ThermoIce0Growth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

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

    class OceanData : public OceanState {
    public:
        OceanData()
            : OceanState()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            sst[0] = -1.75;
            sss[0] = 32.;
            mld[0] = 10.25;
        }
        void updateSpecial(const TimestepTime& tst) override { }
    };
    Module::Module<OceanState>::setExternalImplementation(Module::newImpl<OceanState, OceanData>);

    class AtmosphereData : public AtmosphereState {
    public:
        AtmosphereData()
            : AtmosphereState()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            tair[0] = -12;
            tdew[0] = -12;
            pair[0] = 100000.;
            windSpeed[0] = 5;
            sw_in[0] = 0;
            lw_in[0] = 265;
            snowfall[0] = 1e-3;
        }
        void updateSpecial(const TimestepTime& tst) override { }
    };
    Module::Module<AtmosphereState>::setExternalImplementation(
        Module::newImpl<AtmosphereState, AtmosphereData>);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -9.;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    FiniteElementFluxCalc fefc;
    fefc.configure();
    fefc.setData(ModelState().data);
    fefc.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::Q_OW, RO> qow;
    ModelArrayRef<ModelComponent::SharedArray::Q_IO, RO> qio;
    ModelArrayRef<ModelComponent::SharedArray::Q_IA, RO> qia;
    ModelArrayRef<ModelComponent::SharedArray::DQIA_DT, RO> dqia_dt;
    ModelArrayRef<ModelComponent::SharedArray::SUBLIM, RO> subl;

    double prec = 1e-5;
    REQUIRE(qow[0] == Approx(143.266).epsilon(prec));
    REQUIRE(qio[0] == Approx(73.9465).epsilon(prec));
    REQUIRE(qia[0] == Approx(42.2955).epsilon(prec));
    REQUIRE(dqia_dt[0] == Approx(16.7615).epsilon(prec));
    REQUIRE(subl[0] == Approx(2.15132e-6).epsilon(prec));
}
}
