/*!
 * @file ThermoIce0Temperature_test.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <sstream>

#include "include/IFluxCalculation.hpp"
#include "include/ThermoIce0.hpp"

#include "include/Configurator.hpp"
#include "include/constants.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ThermoIce0Temperature");
TEST_CASE("Melting conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
        {
            registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice[0] = 0.5;
            hice[0] = 0.1 / cice[0]; // Here we are using the true thicknesses
            hsnow[0] = 0.01 / cice[0];
            snow[0] = 0.00;
            tice0[0] = -1;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField snow;
        HField tice0;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState().data);

    UniformOcean oceanData(-1, 32., 4.29151e7 / (Water::rho * Water::cp), 0, 0);
    oceanData.setData(ModelState().data);

    class FluxData : public IFluxCalculation {
    public:
        FluxData()
            : IFluxCalculation()
        {
        }
        std::string getName() const override { return "FluxData"; }

        void setData(const ModelState::DataMap&) override
        {
            qow[0] = -109.923;
            qia[0] = -84.5952;
            dqia_dt[0] = 19.7016;
            subl[0] = -7.3858e-06;
        }

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

        void update(const TimestepTime&) override { }
    } fluxData;
    fluxData.setData(ModelState().data);

    TimestepTime tst = { 0, 600 };
    ThermoIce0 ti0t;
    ti0t.configure();
    ti0t.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::T_ICE, MARBackingStore, RO> tice(
        ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::Q_IC, MARBackingStore, RO> qic(
        ModelComponent::getSharedArray());

    double prec = 1e-5;
    REQUIRE(tice[0] == 0.0);
    REQUIRE(qic[0] == doctest::Approx(-4.60879).epsilon(prec));
}

TEST_CASE("Freezing conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing");
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo");
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
        {
            registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice[0] = 0.5;
            hice[0] = 0.1 / cice[0]; // Here we are using the true thicknesses
            hsnow[0] = 0.01 / cice[0];
            sst[0] = -1.75;
            sss[0] = 32.;
            snow[0] = 1e-3;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;
            tice0[0] = -9.;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity
        HField tice0;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState());

    class FluxData : public IFluxCalculation {
    public:
        FluxData()
            : IFluxCalculation()
        {
        }
        std::string getName() const override { return "FluxData"; }

        void setData(const ModelState::DataMap&) override
        {
            qow[0] = 143.266;
            qio[0] = 73.9465;
            qia[0] = 42.2955;
            dqia_dt[0] = 16.7615;
            subl[0] = 2.15132e-6;
        }

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

        void update(const TimestepTime&) override { }
    } fluxData;

    fluxData.setData(ModelState());

    TimestepTime tst = { 0, 600 };
    ThermoIce0 ti0t;
    ti0t.configure();
    ti0t.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::T_ICE, MARBackingStore, RO> tice(
        ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::Q_IC, MARBackingStore, RO> qic(
        ModelComponent::getSharedArray());

    double prec = 1e-5;
    REQUIRE(tice[0] == doctest::Approx(-8.90443).epsilon(prec));
    REQUIRE(qic[0] == doctest::Approx(44.4839).epsilon(prec));
}
TEST_SUITE_END();

}
