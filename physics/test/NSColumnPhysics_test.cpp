/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include <doctest/doctest.h>
#include <sstream>

#include "include/NSColumnPhysics.hpp"

#include "include/Configurator.hpp"
#include "include/IAtmosphereBoundary.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/NextsimModule.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

extern template class Module::Module<Nextsim::IIceThermodynamics>;
namespace Nextsim {

TEST_SUITE_BEGIN("NSColumnPhysics");
TEST_CASE("New ice formation")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<ILateralIceSpread>::setImplementation("Nextsim::HiblerSpread");
    Module::Module<IIceThermodynamics>::setImplementation("Nextsim::ThermoIce0");
    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    class AtmosphereBoundary : public IAtmosphereBoundary {
    public:
        AtmosphereBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qiaAccessor.getHostRW() = 305.288;
            dqia_dtAccessor.getHostRW() = 4.5036;
            qowAccessor.getHostRW() = 307.546;
            sublAccessor.getHostRW() = 0.; // Seems unlikely…
            penSWAccessor.getHostRW() = 0.;
            snowAccessor.getHostRW() = 0.;
            rainAccessor.getHostRW() = 0.;
            evapAccessor.getHostRW() = 0.; // Seems unlikely…
            uwindAccessor.getHostRW() = 0;
            vwindAccessor.getHostRW() = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW() = 0.5;
            hiceAccessor.getHostRW() = 0.1; // Cell averaged
            hsnowAccessor.getHostRW() = 0; // Cell averaged
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } proData;
    proData.setData(ModelState().data);

    UniformOcean ocnBdy(-1.5, 32., 10.25);
    ocnBdy.setQio(124.689);
    ocnBdy.setData(ModelState().data);

    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor(
        ModelComponent::getStore(), RW, ModelArray::Type::H);
    ModelArrayAccessor<Protected::DAMAGE, RW> oldDamageAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    damageAccessor.getHostRW() = 1;

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1") };
    NSColumnPhysics ig;
    ig.configure();

    HField tsurf;
    tsurf = -2.;
    ModelState::DataMap dataMap = { { tsurfName, tsurf } };
    ig.setData(dataMap);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayAccessor<Shared::NEW_ICE, RO> newiceAccessor(ModelComponent::getStore());
    const HField& newice = newiceAccessor.getHostRO();

    double prec = 1e-5;
    REQUIRE(newice[0] == doctest::Approx(0.0258264).epsilon(prec));
}

TEST_CASE("Melting conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<ILateralIceSpread>::setImplementation("Nextsim::HiblerSpread");
    Module::Module<IIceThermodynamics>::setImplementation("Nextsim::ThermoIce0");

    class AtmosphericBoundary : public IAtmosphereBoundary {
    public:
        AtmosphericBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qiaAccessor.getHostRW() = -84.5952;
            dqia_dtAccessor.getHostRW() = 19.7016;
            qowAccessor.getHostRW() = -109.923;
            sublAccessor.getHostRW() = -7.3858e-06;
            penSWAccessor.getHostRW() = 0.;
            snowAccessor.getHostRW() = 0.;
            rainAccessor.getHostRW() = 0.;
            evapAccessor.getHostRW() = 0.; // Seems unlikely…
            uwindAccessor.getHostRW() = 0;
            vwindAccessor.getHostRW() = 0.;
        }
        std::string getName() const override { return "AtmosphericBoundary"; }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW() = 0.5;
            hiceAccessor.getHostRW() = 0.1; // Cell averaged
            hsnowAccessor.getHostRW() = 0.01; // Cell averaged
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } proData;
    proData.setData(ModelState().data);

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1, 32., 10.25);
    ocnBdy.setQio(53717.8);
    ocnBdy.setData(ModelState().data);

    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor(
        ModelComponent::getStore(), RW, ModelArray::Type::H);
    ModelArrayAccessor<Protected::DAMAGE, RW> oldDamageAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    damageAccessor.getHostRW() = 1;

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    NSColumnPhysics ig;
    ig.configure();

    HField tsurf;
    tsurf = -1.;
    ModelState::DataMap dataMap = { { tsurfName, tsurf } };
    ig.setData(dataMap);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayAccessor<Shared::NEW_ICE, RO> newiceAccessor(ModelComponent::getStore());
    const HField& newice = newiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_ICE_DG, RO> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_SNOW_DG, RO> hsnowAccessor(ModelComponent::getStore());
    const HField& hsnow = hsnowAccessor.getHostRO();

    double prec = 1e-5;
    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE(cice[0] == doctest::Approx(0.368269).epsilon(prec));
    REQUIRE((hice[0]) == doctest::Approx(0.0473078).epsilon(prec));
    REQUIRE((hsnow[0]) == doctest::Approx(0.00720977).epsilon(prec));

    REQUIRE(newice[0] == 0.0);
}

TEST_CASE("Freezing conditions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<ILateralIceSpread>::setImplementation("Nextsim::HiblerSpread");
    Module::Module<IIceThermodynamics>::setImplementation("Nextsim::ThermoIce0");

    class AtmosphereBoundary : public IAtmosphereBoundary {
    public:
        AtmosphereBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qiaAccessor.getHostRW() = 42.2955;
            dqia_dtAccessor.getHostRW() = 16.7615;
            qowAccessor.getHostRW() = 143.266;
            sublAccessor.getHostRW() = 2.15132e-6;
            penSWAccessor.getHostRW() = 0.;
            snowAccessor.getHostRW() = 1e-3;
            rainAccessor.getHostRW() = 0.;
            evapAccessor.getHostRW() = -1e-3; // E-P = 0
            uwindAccessor.getHostRW() = 0;
            vwindAccessor.getHostRW() = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW() = 0.5;
            hiceAccessor.getHostRW() = 0.1; // Cell averaged
            hsnowAccessor.getHostRW() = 0.01; // Cell averaged
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } proData;
    proData.setData(ModelState().data);

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    UniformOcean ocnBdy(-1.75, 32., 10.25);
    ocnBdy.setQio(73.9465);
    ocnBdy.setData(ModelState().data);

    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor(
        ModelComponent::getStore(), RW, ModelArray::Type::H);
    ModelArrayAccessor<Protected::DAMAGE, RW> oldDamageAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    damageAccessor.getHostRW() = 1;

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    NSColumnPhysics ig;
    ig.configure();

    HField tsurf;
    tsurf = -9.;
    ModelState::DataMap dataMap = { { tsurfName, tsurf } };
    ig.setData(dataMap);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayAccessor<Shared::NEW_ICE, RO> newiceAccessor(ModelComponent::getStore());
    const HField& newice = newiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_ICE_DG, RO> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_SNOW_DG, RO> hsnowAccessor(ModelComponent::getStore());
    const HField& hsnow = hsnowAccessor.getHostRO();

    double prec = 1e-5;

    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE(cice[0] == doctest::Approx(0.5002).epsilon(prec));
    REQUIRE((hice[0]) == doctest::Approx(0.100039).epsilon(prec));
    REQUIRE((hsnow[0]) == doctest::Approx(0.0109012).epsilon(prec));

    REQUIRE(newice[0] == doctest::Approx(6.79906e-5).epsilon(prec));
}

TEST_CASE("Dummy ice")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

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
            qiaAccessor.getHostRW() = 0.;
            dqia_dtAccessor.getHostRW() = 0.;
            qowAccessor.getHostRW() = 0.;
            sublAccessor.getHostRW() = 0.;
            penSWAccessor.getHostRW() = 0.;
            snowAccessor.getHostRW() = 0.;
            rainAccessor.getHostRW() = 0.;
            evapAccessor.getHostRW() = 0.; // E-P = 0
            uwindAccessor.getHostRW() = 0;
            vwindAccessor.getHostRW() = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    // Don't like referencing variables in the enclosing scope? FINE!
#define cice0 0.5
#define hice0 0.1
#define hsnow0 0.01
#define tice00 -5

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW() = cice0;
            hiceAccessor.getHostRW() = hice0; // Cell averaged
            hsnowAccessor.getHostRW() = hsnow0; // Cell averaged
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } proData;
    proData.setData(ModelState().data);

    UniformOcean ocnBdy(-1, 35, 10);
    ocnBdy.setQio(0.);
    ocnBdy.setData(ModelState().data);

    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor(
        ModelComponent::getStore(), RW, ModelArray::Type::H);
    ModelArrayAccessor<Protected::DAMAGE, RW> oldDamageAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    damageAccessor.getHostRW() = 1;

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };

    NSColumnPhysics ig;
    ig.configure();

    HField tsurf;
    tsurf = tice00;
    ModelState::DataMap dataMap = { { tsurfName, tsurf } };
    ig.setData(dataMap);
    ocnBdy.updateBefore(tst);

    ig.update(tst);

    //   double prec = 1e-5;

    ModelArrayAccessor<Shared::NEW_ICE, RO> newiceAccessor(ModelComponent::getStore());
    const HField& newice = newiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_ICE_DG, RO> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_SNOW_DG, RO> hsnowAccessor(ModelComponent::getStore());
    const HField& hsnow = hsnowAccessor.getHostRO();

    // The thickness values from old NextSIM are cell-averaged. Perform that
    // conversion here.
    REQUIRE(cice[0] == cice0);
    REQUIRE((hice[0]) == hice0);
    REQUIRE((hsnow[0]) == hsnow0);

    REQUIRE(newice[0] == 0.);
}
#undef cice0
#undef hice0
#undef hsnow0
#undef tice00

TEST_CASE("Zero thickness")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<ILateralIceSpread>::setImplementation("Nextsim::HiblerSpread");
    Module::Module<IIceThermodynamics>::setImplementation("Nextsim::ThermoIce0");

    class AtmosphericBoundary : public IAtmosphereBoundary {
    public:
        AtmosphericBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qiaAccessor.getHostRW() = -84.5952;
            dqia_dtAccessor.getHostRW() = 19.7016;
            qowAccessor.getHostRW() = -109.923;
            sublAccessor.getHostRW() = -7.3858e-06;
            penSWAccessor.getHostRW() = 0.;
            snowAccessor.getHostRW() = 0.;
            rainAccessor.getHostRW() = 0.;
            evapAccessor.getHostRW() = 0.; // Seems unlikely…
            uwindAccessor.getHostRW() = 0;
            vwindAccessor.getHostRW() = 0.;
        }
        std::string getName() const override { return "AtmosphericBoundary"; }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW() = 0.5;
            hiceAccessor.getHostRW() = 0.1; // Cell averaged
            hsnowAccessor.getHostRW() = 0.01; // Cell averaged
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
    } proData;
    proData.setData(ModelState().data);

    UniformOcean ocnBdy(-1, 32., 10.25);
    ocnBdy.setQio(53717.8); // 57 kW m⁻² to go from -1 to -1.75 over the whole mixed layer in 600 s
    ocnBdy.setData(ModelState().data);

    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor(
        ModelComponent::getStore(), RW, ModelArray::Type::H);
    ModelArrayAccessor<Protected::DAMAGE, RW> oldDamageAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    damageAccessor.getHostRW() = 1;

    class ZeroThicknessIce : public IIceThermodynamics {
        void setData(const ModelState::DataMap&) override { }
        void update(const TimestepTime& tsTime) override
        {
            AdvectedField& hice = hiceAccessor.getHostRW();
            deltaHiAccessor.getHostRW()[0] = -hice[0];
            hice[0] = 0;
            tsurfAccessor.getHostRW()[0] = 0;
            snowToIceAccessor.getHostRW()[0] = 0;
        }
    };
    Module::Module<IIceThermodynamics>::setExternalImplementation(
        Module::newImpl<IIceThermodynamics, ZeroThicknessIce>);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1") };
    NSColumnPhysics ig;
    ig.configure();

    HField tsurf;
    tsurf = -1.;
    ModelState::DataMap dataMap = { { tsurfName, tsurf } };
    ig.setData(dataMap);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayAccessor<Shared::NEW_ICE, RO> newiceAccessor(ModelComponent::getStore());
    const HField& newice = newiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_ICE_DG, RO> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();

    //    double prec = 1e-6;

    REQUIRE(newice[0] == 0);
    REQUIRE(hice[0] == 0);
    REQUIRE(cice[0] == 0);
}

TEST_CASE("Turn off thermo")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::Module<ILateralIceSpread>::setImplementation("Nextsim::HiblerSpread");
    Module::Module<IIceThermodynamics>::setImplementation("Nextsim::ThermoIce0");
    std::stringstream config;
    config << "[nextsim_thermo]" << std::endl;
    config << "use_thermo_forcing = false" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    class AtmosphereBoundary : public IAtmosphereBoundary {
    public:
        AtmosphereBoundary()
            : IAtmosphereBoundary()
        {
        }
        void setData(const ModelState::DataMap& ms) override
        {
            IAtmosphereBoundary::setData(ms);
            qiaAccessor.getHostRW() = 42.2955;
            dqia_dtAccessor.getHostRW() = 16.7615;
            qowAccessor.getHostRW() = 143.266;
            sublAccessor.getHostRW() = 2.15132e-6;
            penSWAccessor.getHostRW() = 0.;
            snowAccessor.getHostRW() = 1e-3;
            rainAccessor.getHostRW() = 0.;
            evapAccessor.getHostRW() = -1e-3; // E-P = 0
            uwindAccessor.getHostRW() = 0;
            vwindAccessor.getHostRW() = 0.;
        }
    } atmBdy;
    atmBdy.setData(ModelState().data);

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
            : hiceAccessor(getStore(), RW)
            , ciceAccessor(getStore(), RW)
            , hsnowAccessor(getStore(), RW)
        {
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            ciceAccessor.getHostRW() = 0.5;
            hiceAccessor.getHostRW() = 0.1; // Cell averaged
            hsnowAccessor.getHostRW() = 0.01; // Cell averaged
        }

        ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
        ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
        ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowAccessor;
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
            qioAccessor.getHostRW() = 73.9465;
            sstAccessor.getHostRW() = -1.75;
            sssAccessor.getHostRW() = 32.;
            mldAccessor.getHostRW() = 10.25;
            uAccessor.getHostRW() = 0.;
            vAccessor.getHostRW() = 0.;
        }
        void updateBefore(const TimestepTime& tst) override
        {
            UnescoFreezing uf;
            const HField& mld = mldAccessor.getHostRO();
            cpmlAccessor.getHostRW() = Water::cp * Water::rho * mld;
            const HField& sss = sssAccessor.getHostRO();
            tfAccessor.getHostRW() = uf(sss[0]);
        }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor(
        ModelComponent::getStore(), RW, ModelArray::Type::H);
    ModelArrayAccessor<Protected::DAMAGE, RW> oldDamageAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    damageAccessor.getHostRW() = 1;

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    NSColumnPhysics ig;
    ig.configure();

    HField tsurf;
    tsurf = -9.;
    ModelState::DataMap dataMap = { { tsurfName, tsurf } };
    ig.setData(dataMap);
    ocnBdy.updateBefore(tst);
    ig.update(tst);

    ModelArrayAccessor<Shared::NEW_ICE, RO> newiceAccessor(ModelComponent::getStore());
    const HField& newice = newiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_ICE_DG, RO> hiceAccessor(ModelComponent::getStore());
    const HField& hice = hiceAccessor.getHostRO();
    ModelArrayAccessor<Shared::C_ICE_DG, RO> ciceAccessor(ModelComponent::getStore());
    const HField& cice = ciceAccessor.getHostRO();
    ModelArrayAccessor<Shared::H_SNOW_DG, RO> hsnowAccessor(ModelComponent::getStore());
    const HField& hsnow = hsnowAccessor.getHostRO();

    //    double prec = 1e-5;

    // Rather than the values from old NextSIM, they should be unchanged from the definition above.
    REQUIRE(cice[0] == 0.5);
    REQUIRE((hice[0]) == 0.1);
    REQUIRE((hsnow[0]) == 0.01);

    REQUIRE(newice[0] == 0.0);
}

TEST_SUITE_END();

}
