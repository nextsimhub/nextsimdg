/*!
 * @file ElementData_test.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/ConfiguredModule.hpp"
#include "include/ElementData.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/ModuleLoader.hpp"
#include "include/NextsimPhysics.hpp"

namespace Nextsim {

TEST_CASE("Physics test using NextsimPhysics", "[ElementData]")
{
    // Largely copied from the "Melting conditions" test of NextsimPhysics
    Configurator::clear();
    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;
    config << "Nextsim::IIceAlbedo = Nextsim::CCSMIceAlbedo" << std::endl;
    config << "Nextsim::IPhysics1d = Nextsim::NextsimPhysics" << std::endl;
    config << std::endl;
    config << "[CCSMIceAlbedo]" << std::endl;
    config << "iceAlbedo = 0.63" << std::endl;
    config << "snowAlbedo = 0.88" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    double tair = 3; //˚C
    double tdew = 2; //˚C
    double pair = 100000; // Pa, slightly low pressure
    double sst = -1; //˚C
    double sss = 32; // PSU
    std::vector<double> tice = { -1., -1., -1. }; //˚C
    double hice = 0.1; // m
    double cice = 0.5;
    double hsnow = 0.01; // m
    double dml = 10.; // m

    ModuleLoader::getLoader().setAllDefaults();
    ConfiguredModule::parseConfigurator();
    tryConfigure(ModuleLoader::getLoader().getImplementation<IIceAlbedo>());

    ElementData data(3);
    REQUIRE(data.nIceLayers() == 3);
    data.configure(); // Configure with the UNESCO freezing point

    data =
        PrognosticGenerator().hice(hice).cice(cice).sst(sst).sss(sss).hsnow(hsnow).tice(tice);
    data.setTimestep(600.); // s. Very long TS to get below freezing
    REQUIRE(data.iceThickness() == hice);
    REQUIRE(data.iceConcentration() == cice);
    REQUIRE(data.snowThickness() == hsnow);
    REQUIRE(data.seaSurfaceTemperature() == sst);
    REQUIRE(data.seaSurfaceSalinity() == sss);
    REQUIRE(data.iceTemperature(0) == tice[0]);
    REQUIRE(data.iceTemperature(2) == tice[2]);

    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 330;
    data.incomingShortwave() = 50;
    data.snowfall() = 0;

    data.windSpeed() = 5;

    data.updateDerivedData(data, data, data);
    REQUIRE(1.26488 == Approx(data.airDensity()).epsilon(1e-4));
    REQUIRE(0.00441973 == Approx(data.specificHumidityAir()).epsilon(1e-4));
    REQUIRE(0.0035214 == Approx(data.specificHumidityIce()).epsilon(1e-4));
    REQUIRE(1012.86 == Approx(data.heatCapacityWetAir()).epsilon(1e-4));

    data.calculate(data, data, data);

    REQUIRE(0.12846 == Approx(data.updatedIceTrueThickness()).epsilon(1e-4));
    REQUIRE(0.01957732 == Approx(data.updatedSnowTrueThickness()).epsilon(1e-4));
    REQUIRE(0.368269 == Approx(data.updatedIceConcentration()).epsilon(1e-4));
    REQUIRE(0.0 == Approx(data.updatedIceSurfaceTemperature()).epsilon(1e-4));
}
} /* namespace Nextsim */
