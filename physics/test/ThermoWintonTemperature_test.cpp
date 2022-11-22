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
#include "include/IFreezingPoint.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

#include <iostream> // FIXME remove me

namespace Nextsim {

TEST_CASE("Melting conditions", "[ThermoWinton]")
{
    ModelArray::setDimensions(ModelArray::Type::H,
        {
            1,
        });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 3 });

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

    std::cerr << "configuration done" << std::endl;

    ThermoWinton twin;
    class IceTemperatureData : public ModelComponent {
    public:
        IceTemperatureData()
            : tice0(ModelArray::Type::Z)
            , tice(getSharedArray())
        {
            registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice0);
            registerProtectedArray(ProtectedArray::C_ICE, &cice0);
            registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow0);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::SW_IN, &sw_in);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);

            registerSharedArray(SharedArray::H_ICE, &hice);
            registerSharedArray(SharedArray::C_ICE, &cice);
            registerSharedArray(SharedArray::H_SNOW, &hsnow);
        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice0[0] = 0.5;
            hice0[0] = 0.1 / cice0[0]; // Here we are using the true thicknesses
            hsnow0[0] = 0.01 / cice0[0];
            snow[0] = 0.00;
            sw_in[0] = -50;
            sst[0] = -1;
            sss[0] = 32.;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;
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
        HField snow;
        HField sw_in;
        HField sst;
        HField sss;
        HField tf;
        HField mlbhc; // Mixed layer bulk heat capacity
        ZField tice0;
        ModelArrayRef<SharedArray::T_ICE, MARBackingStore, RW> tice; // From IIceThermodynamics

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState().data);

    std::cerr << "atmos data defined" << std::endl;
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
            qio[0]
                = 53717.8; // 57 kW m⁻² to go from -1 to -1.75 over the whole mixed layer in 600 s
            qia[0] = -84.5952;
            dqia_dt[0] = 19.7016;
            subl[0] = -7.3858e-06;
        }

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

        void update(const TimestepTime&) override { }
    } fluxData;
    fluxData.setData(ModelState().data);

    std::cerr << "flux data defined" << std::endl;
    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    twin.configure();
    twin.update(tst);
    std::cerr << "thermo winton updated" << std::endl;
    ModelArrayRef<ModelComponent::SharedArray::T_ICE, MARBackingStore, RO> tice(
        ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::Q_IC, MARBackingStore, RO> qic(
        ModelComponent::getSharedArray());

    double prec = 1e-5;
    REQUIRE(tice[0] == Approx(0.0).margin(prec));
    REQUIRE(tice[1] == Approx(0.0).margin(prec));
    REQUIRE(tice[2] == Approx(-0.999453).epsilon(prec));
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
        , tice(getSharedArray())
        {
            registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice0);
            registerProtectedArray(ProtectedArray::C_ICE, &cice0);
            registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow0);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::SW_IN, &sw_in);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);

            registerSharedArray(SharedArray::H_ICE, &hice);
            registerSharedArray(SharedArray::C_ICE, &cice);
            registerSharedArray(SharedArray::H_SNOW, &hsnow);

        }
        std::string getName() const override { return "IceTemperatureData"; }

        void setData(const ModelState::DataMap&) override
        {
            cice0[0] = 0.5;
            hice0[0] = 0.1 / cice0[0]; // Here we are using the true thicknesses
            hsnow0[0] = 0.01 / cice0[0];
            snow[0] = 1e-3;
            sw_in[0] = -50;
            sst[0] = -1.75;
            sss[0] = 32.;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;
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
        HField sst;
        HField sss;
        HField tf;
        HField mlbhc; // Mixed layer bulk heat capacity
        ZField tice0;
        ModelArrayRef<SharedArray::T_ICE, MARBackingStore, RW> tice; // From IIceThermodynamics

        HField hice;
        HField cice;
        HField hsnow;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } atmoState;
    atmoState.setData(ModelState().data);

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

            std::cerr << "FluxData::setData &subl=" << &subl;
            ModelArrayRef<SharedArray::SUBLIM, MARBackingStore, RO> sublRef(getSharedArray());
            std::cerr << " &SharedArray::SUBLIM.data()=" << &sublRef.data() << std::endl;
}

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

        void update(const TimestepTime&) override { }
    } fluxData;

    fluxData.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    twin.configure();
    ModelArrayRef<ModelComponent::SharedArray::SUBLIM, MARBackingStore, RO> sublRef(ModelComponent::getSharedArray());
    std::cerr << " &SharedArray::SUBLIM.data()=" << &sublRef.data() << std::endl;

    twin.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::T_ICE, MARBackingStore, RO> tice(
        ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::Q_IC, MARBackingStore, RO> qic(
        ModelComponent::getSharedArray());

    double prec = 1e-5;
    REQUIRE(tice[0] == Approx(-8.90443).epsilon(prec));
    REQUIRE(qic[0] == Approx(44.4839).epsilon(prec));
}
}
