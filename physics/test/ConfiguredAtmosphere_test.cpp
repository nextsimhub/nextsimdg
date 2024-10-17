/*!
 * @file ConfiguredAtmosphere_test.cpp
 *
 * @date Sep 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <sstream>

#include "include/ConfiguredAtmosphere.hpp"

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/IFluxCalculation.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/NextsimModule.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/UniformOcean.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ConfiguredAtmosphere");
TEST_CASE("ConfiguredAtmosphere melting test")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    Module::Module<IFreezingPoint>::setImplementation("Nextsim::UnescoFreezing")
    Module::Module<IIceAlbedo>::setImplementation("Nextsim::CCSMIceAlbedo")
    Module::Module<IFluxCalculation>::setImplementation("Nextsim::FiniteElementFluxes")
    std::stringstream config;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;
    config << "[ConfiguredAtmosphere]" << std::endl;
    config << "t_air = 3" << std::endl;
    config << "t_dew = 2" << std::endl;
    config << "pmsl = 100000" << std::endl;
    config << "sw_in = 50" << std::endl;
    config << "lw_in = 330" << std::endl;
    config << "snow = 0" << std::endl;
    config << "rainfall = 0" << std::endl;
    config << "wind_speed = 5" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");
    UniformOcean ocnBdy(-1., 32., 10.25);
    ocnBdy.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            getStore().registerArray(Protected::H_ICE, &hice, RO);
            getStore().registerArray(Protected::C_ICE, &cice, RO);
            getStore().registerArray(Protected::H_SNOW, &hsnow, RO);
            getStore().registerArray(Protected::T_ICE, &tice0, RO);
            getStore().registerArray(Protected::HTRUE_ICE, &hice0, RO);
            getStore().registerArray(Protected::HTRUE_SNOW, &hsnow0, RO);
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -1.;

            hice0[0] = hice[0] / cice[0];
            hsnow0[0] = hsnow[0] / cice[0];
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        HField hice0; // ice averaged ice thickness
        HField hsnow0; // ice averaged snow thickness
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    ConfiguredAtmosphere ca;
    ca.configure();
    ca.setData(ModelState().data);

    HField qow;
    qow.resize();
    ModelComponent::getStore().registerArray(Shared::Q_OW, &qow, RW);

    HField qia;
    qia.resize();
    ModelComponent::getStore().registerArray(Shared::Q_IA, &qia, RW);

    HField dqia_dt;
    dqia_dt.resize();
    ModelComponent::getStore().registerArray(Shared::DQIA_DT, &dqia_dt, RW);

    HField subl;
    subl.resize();
    ModelComponent::getStore().registerArray(Shared::SUBLIM, &subl, RW);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    ocnBdy.updateBefore(tst);
    ca.update(tst);

    double prec = 1e-5;
    REQUIRE(qow[0] == doctest::Approx(-109.923).epsilon(prec));
    REQUIRE(qia[0] == doctest::Approx(-85.6364).epsilon(prec));
    REQUIRE(dqia_dt[0] == doctest::Approx(19.7016).epsilon(prec));
    REQUIRE(subl[0] == doctest::Approx(-7.3858e-06).epsilon(prec));
}
TEST_SUITE_END();
}
