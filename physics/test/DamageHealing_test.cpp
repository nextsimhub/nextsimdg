/*!
 * @file DamageHealing_test.cpp
 *
 * @date Jul 4, 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/NextsimModule.hpp"
#include "include/IDamageHealing.hpp"

extern template class Module::Module<Nextsim::IDamageHealing>;
namespace Nextsim {

TEST_SUITE_BEGIN("DamageHealing");
TEST_CASE("Thermodynamic healing")
/* Test damage healing due to thermodynamics
 * We set the healing timescale (td) to a sensible value and impose no new-ice growth.
 * We expect damage to grow linearly, proportional to dt/td, where dt is the time step.
 *
 * The tests here are:
 * damage = 0.5 & td = 20 days & dt = 1 day => damage = 0.55
 * damage = 0.99 & td = 20 days & dt = 1 day => damage = 1.0 (must not exceed 1)
 */
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    Module::Module<IDamageHealing>::setImplementation("Nextsim::ConstantHealing");
    std::stringstream config;

    config << "[ConstantHealing]" << std::endl;
    config << "td = 20" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
        {
            getStore().registerArray(Shared::DELTA_CICE, &deltaCi, RO);
            getStore().registerArray(Shared::C_ICE, &cice, RO);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            deltaCi= 0.0;
            damage = 0.5;
        }

        HField cice;
        HField deltaCi;
        HField damage;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    std::unique_ptr<IDamageHealing> iHealing;
    iHealing = std::move(Module::getInstance<IDamageHealing>());
    tryConfigure(*iHealing);
    iHealing->setData(ModelState().data);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1T00:00:00") };
    double prec = 1e-8;

    iceState.damage = 0.5;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] == doctest::Approx(0.55).epsilon(prec));

    iceState.damage = 0.99;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] <= 1. );
    REQUIRE(iceState.damage[0] == doctest::Approx(1.).epsilon(prec));
}

TEST_CASE("New ice formation")
/* Test damage healing due to new ice formation.
 * We expect damage to grow proportionally to the amount of ice area added, i.e. all new, laterally
 * grown ice has damage = 1.
 *
 * cice is the ice concentration after adding deltaCi damage is the damage
 * before any healing takes place
 *
 * The tests here are
 * cice = 0.6 & deltaCi = 0.3 & damage 0 => damage = 0.5
 * cice = 0.6 & deltaCi = 0.3 & damage 0.5 => damage = 0.75
 * cice = 0.6 & deltaCi = 0.1 & damage 0.5 => damage = 0.6
 * cice = 1.0 & deltaCi = 0.1 & damage 1.0 => damage = 1.0 (must not exceed 1)
 * Melting (deltaCi < 0) should not change the level of damage
 *
 * I must always add 0.05, because I can't reconfigure the module to change the healing timescale
 * (see previous test).
 */
{
    class PrognosticData : public ModelComponent {
    public:
        PrognosticData()
        {
            getStore().registerArray(Shared::DELTA_CICE, &deltaCi, RO);
            getStore().registerArray(Shared::C_ICE, &cice, RO);
            getStore().registerArray(Shared::DAMAGE, &damage, RW);
        }
        std::string getName() const override { return "PrognosticData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice = 0.5;
            deltaCi= 0.1;
            damage = 0.5;
        }

        HField cice;
        HField deltaCi;
        HField damage;

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    std::unique_ptr<IDamageHealing> iHealing;
    iHealing = std::move(Module::getInstance<IDamageHealing>());

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-1T00:00:00") };
    double prec = 1e-8;

    iceState.cice = 0.6;
    iceState.deltaCi = 0.3;
    iceState.damage = 0.;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] == doctest::Approx(0.55).epsilon(prec));

    iceState.cice = 0.6;
    iceState.deltaCi = 0.3;
    iceState.damage = 0.5;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] == doctest::Approx(0.80).epsilon(prec));

    iceState.cice = 0.5;
    iceState.deltaCi = 0.1;
    iceState.damage = 0.5;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] == doctest::Approx(0.65).epsilon(prec));

    iceState.cice = 1.;
    iceState.deltaCi = 0.1;
    iceState.damage = 1.;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] <= 1.);
    REQUIRE(iceState.damage[0] <= doctest::Approx(1.).epsilon(prec));

    iceState.cice = 0.5;
    iceState.deltaCi = -0.5;
    iceState.damage = 0.5;
    iHealing->update(tst);
    REQUIRE(iceState.damage[0] == doctest::Approx(0.55).epsilon(prec));
}
}
