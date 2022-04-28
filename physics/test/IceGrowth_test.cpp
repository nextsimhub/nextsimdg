/*!
 * @file IceGrowth_test.cpp
 *
 * @date Apr 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/IceGrowth.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

TEST_CASE("New ice formation", "[IceGrowth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class NewIceData : public ModelComponent {
    public:
        NewIceData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);

            registerSharedArray(SharedArray::Q_OW, &qow);
            registerSharedArray(SharedArray::Q_IC, &qic);
            registerSharedArray(SharedArray::Q_IO, &qio);
            registerSharedArray(SharedArray::Q_IA, &qia);
            registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
            registerSharedArray(SharedArray::SUBLIM, &subl);

            registerSharedArray(SharedArray::T_ICE, &tice);
        }
        std::string getName() const override { return "NewIceData"; }

        void setData(const ModelState&) override
        {
            hice[0] = 0.2;
            cice[0] = 0.5;
            hsnow[0] = 0;
            sst[0] = -1.5;
            sss[0] = 32.;
            snow[0] = 0.;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;

            qow[0] = 307.546;
            qio[0] = 124.689;
            qic[0] = 2.53124;
            qia[0] = 305.288;
            dqia_dt[0] = 4.5036;
            subl[0] = 0;

            tice[0] = -2;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity

        HField qow;
        HField qic;
        HField qio;
        HField qia;
        HField dqia_dt;
        HField subl;

        ZField tice;
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

    } iceData;

    iceData.setData(ModelState());

    TimestepTime tst = {0, 86400};
    IceGrowth ig;
    ig.configure();
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, RO> newice;

    double prec = 1e-5;
    REQUIRE(newice[0] == Approx(0.0258264).epsilon(prec));
}

TEST_CASE("Melting conditions", "[IceGrowth]")
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

    class NewIceData : public ModelComponent {
    public:
        NewIceData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);

            registerSharedArray(SharedArray::Q_OW, &qow);
            registerSharedArray(SharedArray::Q_IC, &qic);
            registerSharedArray(SharedArray::Q_IO, &qio);
            registerSharedArray(SharedArray::Q_IA, &qia);
            registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
            registerSharedArray(SharedArray::SUBLIM, &subl);

            registerSharedArray(SharedArray::T_ICE, &tice);
        }
        std::string getName() const override { return "NewIceData"; }

        void setData(const ModelState&) override
        {
            cice[0] = 0.5;
            hice[0] = 0.1;
            hsnow[0] = 0.01;
            sst[0] = -1;
            sss[0] = 32.;
            snow[0] = 0.00;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;

            qow[0] = -109.923;
            qio[0] = 53717.8; // 57 kW m⁻² to go from -1 to -1.75 over the whole mixed layer in 600 s
            qic[0] = -4.60879;
            qia[0] = -84.5952;
            dqia_dt[0] = 19.7016;
            subl[0] = -7.3858e-06;

            tice[0] = -1;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity

        HField qow;
        HField qic;
        HField qio;
        HField qia;
        HField dqia_dt;
        HField subl;

        ZField tice;
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

    } iceData;

    iceData.setData(ModelState());

    TimestepTime tst = {0, 600};
    IceGrowth ig;
    ig.configure();
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, RO> newice;
    ModelArrayRef<ModelComponent::SharedArray::H_ICE, RO> hice;
    ModelArrayRef<ModelComponent::SharedArray::C_ICE, RO> cice;
    ModelArrayRef<ModelComponent::SharedArray::H_SNOW, RO> hsnow;

    double prec = 1e-5;
    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE((hice[0] * cice[0]) == Approx(0.0473078).epsilon(prec));
    REQUIRE((hsnow[0] * cice[0]) == Approx(0.00720977).epsilon(prec));
    REQUIRE(cice[0] == Approx(0.368269).epsilon(prec));

    REQUIRE(newice[0] == 0.0);
}

TEST_CASE("Freezing conditions", "[IceGrowth]")
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

    class NewIceData : public ModelComponent {
    public:
        NewIceData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::SST, &sst);
            registerProtectedArray(ProtectedArray::SSS, &sss);
            registerProtectedArray(ProtectedArray::TF, &tf);
            registerProtectedArray(ProtectedArray::SNOW, &snow);
            registerProtectedArray(ProtectedArray::ML_BULK_CP, &mlbhc);

            registerSharedArray(SharedArray::Q_OW, &qow);
            registerSharedArray(SharedArray::Q_IC, &qic);
            registerSharedArray(SharedArray::Q_IO, &qio);
            registerSharedArray(SharedArray::Q_IA, &qia);
            registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
            registerSharedArray(SharedArray::SUBLIM, &subl);

            registerSharedArray(SharedArray::T_ICE, &tice);
        }
        std::string getName() const override { return "NewIceData"; }

        void setData(const ModelState&) override
        {
            cice[0] = 0.5;
            hice[0] = 0.1;
            hsnow[0] = 0.01;
            sst[0] = -1.75;
            sss[0] = 32.;
            snow[0] = 1e-3;
            tf[0] = Module::getImplementation<IFreezingPoint>()(sss[0]);
            mlbhc[0] = 4.29151e7;

            qow[0] = 143.266;
            qio[0] = 73.9465;
            qic[0] = 44.4839;
            qia[0] = 42.2955;
            dqia_dt[0] = 16.7615;
            subl[0] = 2.15132e-6;

            tice[0] = -9;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField sst;
        HField sss;
        HField tf;
        HField snow;
        HField mlbhc; // Mixed layer bulk heat capacity

        HField qow;
        HField qic;
        HField qio;
        HField qia;
        HField dqia_dt;
        HField subl;

        ZField tice;
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }

    } iceData;

    iceData.setData(ModelState());

    TimestepTime tst = {0, 600};
    IceGrowth ig;
    ig.configure();
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, RO> newice;
    ModelArrayRef<ModelComponent::SharedArray::H_ICE, RO> hice;
    ModelArrayRef<ModelComponent::SharedArray::C_ICE, RO> cice;
    ModelArrayRef<ModelComponent::SharedArray::H_SNOW, RO> hsnow;

    double prec = 1e-5;

    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE((hice[0] * cice[0]) == Approx(0.100039).epsilon(prec));
    REQUIRE((hsnow[0] * cice[0]) == Approx(0.0109012).epsilon(prec));
    REQUIRE(cice[0] == Approx(0.5002).epsilon(prec));

    REQUIRE(newice[0] == Approx(6.79906e-5).epsilon(prec));
}
}
