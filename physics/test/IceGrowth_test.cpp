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
#include "include/IAtmosphereBoundary.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/IOceanBoundary.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_CASE("New ice formation", "[IceGrowth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::ILateralIceSpread = Nextsim::HiblerSpread" << std::endl;
    config << "Nextsim::IIceThermodynamics = Nextsim::ThermoIce0" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class AtmosphereBoundary : public IAtmosphereBoundary {
    public:
        AtmosphereBoundary()
        : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qia = 305.288;
            dqia_dt = 4.5036;
            qow = 307.546;
            subl = 0.; // Seems unlikely…
            snow = 0.;
            rain = 0.;
            evap = 0.; // Seems unlikely…
            uwind = 0;
            vwind = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0; // Cell averaged
            tice0 = -2;
        }

        HField hice;
        HField cice;
        HField hsnow;
        ZField tice0;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    class OceanBoundary : public IOceanBoundary {
    public:
        OceanBoundary()
            : IOceanBoundary()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            IOceanBoundary::setData(state);
            qio = 124.689;
            sst = -1.5;
            sss = 32.;
            mld = 10.25;
            u = 0.;
            v = 0.;
        }
        void updateBefore(const TimestepTime& tst) override
        {
            UnescoFreezing uf;
            cpml = Water::cp * Water::rho * mld;
            tf = uf(sss[0]);
        }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, MARBackingStore, RO> newice(ModelComponent::getSharedArray());

    double prec = 1e-5;
    REQUIRE(newice[0] == Approx(0.0258264).epsilon(prec));
}

TEST_CASE("Melting conditions", "[IceGrowth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::ILateralIceSpread = Nextsim::HiblerSpread" << std::endl;
    config << "Nextsim::IIceThermodynamics = Nextsim::ThermoIce0" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class AtmosphericBoundary : public IAtmosphereBoundary {
    public:
        AtmosphericBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qia = -84.5952;
            dqia_dt = 19.7016;
            qow = -109.923;
            subl = -7.3858e-06;
            snow = 0.;
            rain = 0.;
            evap = 0.; // Seems unlikely…
            uwind = 0;
            vwind = 0.;
        }
        std::string getName() const override { return "AtmosphericBoundary"; }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0.01; // Cell averaged
            tice0 = -1;
        }

        HField hice;
        HField cice;
        HField hsnow;
        ZField tice0;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    class OceanBoundary : public IOceanBoundary {
    public:
        OceanBoundary()
            : IOceanBoundary()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            IOceanBoundary::setData(state);
            qio = 53717.8; // 57 kW m⁻² to go from -1 to -1.75 over the whole mixed layer in 600 s
            sst[0] = -1;
            sss[0] = 32.;
            mld[0] = 10.25;
            u = 0;
            v = 0;
        }
        void updateBefore(const TimestepTime& tst) override
        {
            UnescoFreezing uf;
            cpml = Water::cp * Water::rho * mld;
            tf = uf(sss[0]);
        }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, MARBackingStore, RO> newice(ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::H_ICE, MARBackingStore, RO> hice(ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::C_ICE, MARBackingStore, RO> cice(ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::H_SNOW, MARBackingStore, RO> hsnow(ModelComponent::getSharedArray());

    double prec = 1e-5;
    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE(cice[0] == Approx(0.368269).epsilon(prec));
    REQUIRE((hice[0] * cice[0]) == Approx(0.0473078).epsilon(prec));
    REQUIRE((hsnow[0] * cice[0]) == Approx(0.00720977).epsilon(prec));

    REQUIRE(newice[0] == 0.0);
}

TEST_CASE("Freezing conditions", "[IceGrowth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::ILateralIceSpread = Nextsim::HiblerSpread" << std::endl;
    config << "Nextsim::IIceThermodynamics = Nextsim::ThermoIce0" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    class AtmosphereBoundary : public IAtmosphereBoundary {
    public:
        AtmosphereBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qia = 42.2955;
            dqia_dt = 16.7615;
            qow = 143.266;
            subl = 2.15132e-6;
            snow = 1e-3;
            rain = 0.;
            evap = -1e-3; // E-P = 0
            uwind = 0;
            vwind = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0.01; // Cell averaged
            tice0 = -9;
        }

        HField hice;
        HField cice;
        HField hsnow;
        ZField tice0;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    class OceanBoundary : public IOceanBoundary {
    public:
        OceanBoundary()
            : IOceanBoundary()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            qio = 73.9465;
            sst = -1.75;
            sss = 32.;
            mld = 10.25;
            u = 0.;
            v = 0.;
        }
        void updateBefore(const TimestepTime& tst) override
        {
            UnescoFreezing uf;
            cpml = Water::cp * Water::rho * mld;
            tf = uf(sss[0]);
        }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<ModelComponent::SharedArray::NEW_ICE, MARBackingStore, RO> newice(ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::H_ICE, MARBackingStore, RO> hice(ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::C_ICE, MARBackingStore, RO> cice(ModelComponent::getSharedArray());
    ModelArrayRef<ModelComponent::SharedArray::H_SNOW, MARBackingStore, RO> hsnow(ModelComponent::getSharedArray());

    double prec = 1e-5;

    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE(cice[0] == Approx(0.5002).epsilon(prec));
    REQUIRE((hice[0] * cice[0]) == Approx(0.100039).epsilon(prec));
    REQUIRE((hsnow[0] * cice[0]) == Approx(0.0109012).epsilon(prec));

    REQUIRE(newice[0] == Approx(6.79906e-5).epsilon(prec));
}
}
