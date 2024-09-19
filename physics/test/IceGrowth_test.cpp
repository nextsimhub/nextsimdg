/*!
 * @file IceGrowth_test.cpp
 *
 * @date 19 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <sstream>

#include "include/IceGrowth.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IAtmosphereBoundary.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Module.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

extern template class Module::Module<Nextsim::IIceThermodynamics>;
namespace Nextsim {

TEST_SUITE_BEGIN("IceGrowth");
TEST_CASE("New ice formation")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "LateralIceSpreadModule = Nextsim::HiblerSpread" << std::endl;
    config << "IceThermodynamicsModule = Nextsim::ThermoIce0" << std::endl;

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
            : hice(ModelArray::Type::H)
            , cice(ModelArray::Type::H)
            , hsnow(ModelArray::Type::H)
            , damage(ModelArray::Type::H)
            , tice(ModelArray::Type::Z)
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::T_ICE, &tice, RW);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0; // Cell averaged
            damage = 1.;
            tice = -2;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField damage;
        ZField tice;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1.5, 32., 10.25);
    ocnBdy.setQio(124.689);
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<Shared::NEW_ICE, RO> newice(ModelComponent::getStore());

    double prec = 1e-5;
    REQUIRE(newice[0] == doctest::Approx(0.0258264).epsilon(prec));
}

TEST_CASE("Melting conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "LateralIceSpreadModule = Nextsim::HiblerSpread" << std::endl;
    config << "IceThermodynamicsModule = Nextsim::ThermoIce0" << std::endl;

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
            : hice(ModelArray::Type::H)
            , cice(ModelArray::Type::H)
            , hsnow(ModelArray::Type::H)
            , damage(ModelArray::Type::H)
            , tice(ModelArray::Type::Z)
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::T_ICE, &tice, RW);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0.01; // Cell averaged
            damage = 1.;
            tice = -1;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField damage;
        ZField tice;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1, 32., 10.25);
    ocnBdy.setQio(53717.8);
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<Shared::NEW_ICE, RW> newice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_ICE, RW> hice(ModelComponent::getStore());
    ModelArrayRef<Shared::C_ICE, RW> cice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_SNOW, RW> hsnow(ModelComponent::getStore());

    double prec = 1e-5;
    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    CHECK(cice[0] == doctest::Approx(0.368269).epsilon(prec));
    CHECK(hice[0] == doctest::Approx(0.0473078).epsilon(prec));
    CHECK(hsnow[0] == doctest::Approx(0.00720977).epsilon(prec));

    CHECK(newice[0] == 0.0);
}

TEST_CASE("Freezing conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "LateralIceSpreadModule = Nextsim::HiblerSpread" << std::endl;
    config << "IceThermodynamicsModule = Nextsim::ThermoIce0" << std::endl;

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
            : hice(ModelArray::Type::H)
            , cice(ModelArray::Type::H)
            , hsnow(ModelArray::Type::H)
            , damage(ModelArray::Type::H)
            , tice(ModelArray::Type::Z)
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::T_ICE, &tice, RW);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0.01; // Cell averaged
            damage = 1.;
            tice = -9;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField damage;
        ZField tice;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1.75, 32., 10.25);
    ocnBdy.setQio(73.9465);
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<Shared::NEW_ICE, RO> newice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_ICE, RO> hice(ModelComponent::getStore());
    ModelArrayRef<Shared::C_ICE, RO> cice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_SNOW, RO> hsnow(ModelComponent::getStore());

    double prec = 1e-5;

    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    CHECK(cice[0] == doctest::Approx(0.5002).epsilon(prec));
    CHECK(hice[0] == doctest::Approx(0.100039).epsilon(prec));
    CHECK(hsnow[0] == doctest::Approx(0.0109012).epsilon(prec));

    CHECK(newice[0] == doctest::Approx(6.79906e-5).epsilon(prec));
}

TEST_CASE("Dummy ice")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    Module::setImplementation<ILateralIceSpread>("Nextsim::DummyIceSpread");
    Module::setImplementation<IIceThermodynamics>("Nextsim::DummyIceThermodynamics");

    class AtmosphereBoundary : public IAtmosphereBoundary {
    public:
        AtmosphereBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qia = 0.;
            dqia_dt = 0.;
            qow = 0.;
            subl = 0.;
            snow = 0.;
            rain = 0.;
            evap = 0.; // E-P = 0
            uwind = 0;
            vwind = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    // Don't like referencing variables in the enclosing scope? FINE!
#define cice0 0.5
#define hice0 0.1
#define hsnow0 0.01
#define damage0 1.
#define tice0 -5

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hice(ModelArray::Type::H)
            , cice(ModelArray::Type::H)
            , hsnow(ModelArray::Type::H)
            , damage(ModelArray::Type::H)
            , tice(ModelArray::Type::Z)
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::T_ICE, &tice, RW);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = cice0;
            hice = hice0; // Cell averaged
            hsnow = hsnow0; // Cell averaged
            damage = damage0;
            tice = tice0;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField damage;
        ZField tice;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    UniformOcean ocnBdy(-1, 35, 10);
    ocnBdy.setQio(0.);
    ocnBdy.setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);

    ig.update(tst);

    double prec = 1e-5;

    ModelArrayRef<Shared::NEW_ICE, RO> newice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_ICE, RO> hice(ModelComponent::getStore());
    ModelArrayRef<Shared::C_ICE, RO> cice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_SNOW, RO> hsnow(ModelComponent::getStore());

    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE(cice[0] == cice0);
    REQUIRE(hice[0] == hice0);
    REQUIRE(hsnow[0] == hsnow0);

    REQUIRE(newice[0] == 0.);
}
#undef cice0
#undef hice0
#undef hsnow0
#undef tice00

TEST_CASE("Zero thickness")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "LateralIceSpreadModule = Nextsim::HiblerSpread" << std::endl;
    config << "IceThermodynamicsModule = Nextsim::ThermoIce0" << std::endl;

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
            : hice(ModelArray::Type::H)
            , cice(ModelArray::Type::H)
            , hsnow(ModelArray::Type::H)
            , damage(ModelArray::Type::H)
            , tice(ModelArray::Type::Z)
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::T_ICE, &tice, RW);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0.01; // Cell averaged
            damage = 1;
            tice = -1;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField damage;
        ZField tice;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } proData;
    proData.setData(ModelState().data);

    UniformOcean ocnBdy(-1, 32., 10.25);
    ocnBdy.setQio(53717.8); // 57 kW m⁻² to go from -1 to -1.75 over the whole mixed layer in 600 s
    ocnBdy.setData(ModelState().data);

    class ZeroThicknessIce : public IIceThermodynamics {
        void setData(const ModelState::DataMap&) override { }
        void update(const TimestepTime& tsTime) override
        {
            deltaHi[0] = -hice[0];
            hice[0] = 0;
            tice[0] = 0;
            snowToIce[0] = 0;
        }
        size_t getNZLevels() const override { return 1; }
    };
    Module::Module<IIceThermodynamics>::setExternalImplementation(
        Module::newImpl<IIceThermodynamics, ZeroThicknessIce>);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1") };
    IceGrowth ig;
    ig.configure();
    ig.setData(ModelState().data);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayRef<Shared::NEW_ICE, RO> newice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_ICE, RO> hice(ModelComponent::getStore());
    ModelArrayRef<Shared::C_ICE, RO> cice(ModelComponent::getStore());

    double prec = 1e-6;

    REQUIRE(newice[0] == 0);
    REQUIRE(hice[0] == 0);
    REQUIRE(cice[0] == 0);
}

TEST_CASE("Turn off thermo")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "LateralIceSpreadModule = Nextsim::HiblerSpread" << std::endl;
    config << "IceThermodynamicsModule = Nextsim::ThermoIce0" << std::endl;
    config << std::endl;
    config << "[nextsim_thermo]" << std::endl;
    config << "use_thermo_forcing = false" << std::endl;

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
            : hice(ModelArray::Type::H)
            , cice(ModelArray::Type::H)
            , hsnow(ModelArray::Type::H)
            , damage(ModelArray::Type::H)
            , tice(ModelArray::Type::Z)
        {
            getStore().registerArray(Shared::H_ICE, &hice, RW);
            getStore().registerArray(Shared::C_ICE, &cice, RW);
            getStore().registerArray(Shared::H_SNOW, &hsnow, RW);
            getStore().registerArray(Shared::T_ICE, &tice, RW);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            hice = 0.1; // Cell averaged
            hsnow = 0.01; // Cell averaged
            damage = 1.;
            tice = -9;
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField damage;
        ZField tice;

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

    ModelArrayRef<Shared::NEW_ICE, RW> newice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_ICE, RW> hice(ModelComponent::getStore());
    ModelArrayRef<Shared::C_ICE, RW> cice(ModelComponent::getStore());
    ModelArrayRef<Shared::H_SNOW, RW> hsnow(ModelComponent::getStore());

    double prec = 1e-5;

    // Rather than the values from old NextSIM, they should be unchanged from the definition above.
    REQUIRE(cice[0] == 0.5);
    REQUIRE(hice[0] == 0.1);
    REQUIRE(hsnow[0] == 0.01);

    REQUIRE(newice[0] == 0.0);
}

TEST_SUITE_END();

}
