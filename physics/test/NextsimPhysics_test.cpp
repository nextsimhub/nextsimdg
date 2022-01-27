/*!
 * @file NextsimPhysics_test.cpp
 *
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
#include "include/constants.hpp"

namespace Nextsim {

TEST_CASE("Minimum ice & i0", "[NextsimPhysics]")
{
    Configurator::clearStreams();

    // Non-default values for the minimum concentration and thickness and I0
    const double minConc = 2e-12;
    const double minThck = 0.02;
    const double i0 = 0.18;

    std::stringstream config;
    config << "[nextsim_thermo]" << std::endl;
    config << "min_conc = " << minConc << std::endl;
    config << "min_thick = " << minThck << std::endl;
    config << "I_0 = " << i0 << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    NextsimPhysics nsphys;
    nsphys.configure();

    REQUIRE(NextsimPhysics::minimumIceConcentration() == minConc);
    REQUIRE(NextsimPhysics::minimumIceThickness() == minThck);
    REQUIRE(NextsimPhysics::i0() == i0);
}

TEST_CASE("Update derived data", "[NextsimPhysics]")
{
    ModuleLoader::getLoader().setAllDefaults();
    ConfiguredModule::parseConfigurator();

    ElementData data;
    data.configure(); // Configure with the default linear freezing point


    double tair = -3;
    double tdew = 0.1;
    double pair = 100000; // Slightly low pressure
    double sst = -1;
    double sss = 32; // PSU
    std::vector<double> tice = { -2., -2, -2 };
    double hice = 0.1;
    double cice = 0.5;

    data = PrognosticData::generate(hice, cice, sst, sss, 0., tice);
    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;

    NextsimPhysics nsData;
    nsData.updateDerivedData(data, data, data);

    REQUIRE(1.29253 == Approx(data.airDensity()).epsilon(1e-4));
    REQUIRE(0.00385326 == Approx(data.specificHumidityAir()).epsilon(1e-4));
    REQUIRE(0.00349446 == Approx(data.specificHumidityWater()).epsilon(1e-4));
    REQUIRE(0.00323958 == Approx(data.specificHumidityIce()).epsilon(1e-4));
    REQUIRE(1011.81 == Approx(data.heatCapacityWetAir()).epsilon(1e-4));
}

TEST_CASE("New ice formation", "[NextsimPhysics]")
{
    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    double tair = -3; //˚C
    double tdew = 0.1; //˚C
    double pair = 100000; // Pa, slightly low pressure
    double sst = -1.5; //˚C
    double sss = 32; // PSU
    std::vector<double> tice = { -2., -2, -2 }; //˚C
    double hice = 0.1; // m
    double cice = 0.5;
    double dml = 10.; // m

    ModuleLoader::getLoader().setAllDefaults();
    ConfiguredModule::parseConfigurator();

    ElementData data(3);
    data.configure(); // Configure with the default linear freezing point

    data = PrognosticData::generate(hice, cice, sst, sss, 0., tice);
    data.setTimestep(86400.); // s. Very long TS to get below freezing

    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 0;
    data.incomingShortwave() = 0;

    NextsimPhysics nsphys;
    nsphys.configure();
    nsphys.updateDerivedData(data, data, data);
    nsphys.calculate(data, data, data);

    // Correct for the non-NIST value of the Stefan-Boltzman constant used by old NeXtSIM. This
    // propagates linearly through to the new sea ice value in this case.
    double sbCorr = PhysicalConstants::sigma / 5.67e-8;
    REQUIRE(0.0258236 * sbCorr == Approx(nsphys.newIce()).epsilon(1e-4));
}

TEST_CASE("Drag pressure", "[NextsimPhysics]")
{
    double tair = 2; //˚C
    double tdew = 1.5; //˚C
    double pair = 100000; // Pa, slightly low pressure
    double sst = -1.5; //˚C
    double sss = 32; // PSU
    std::vector<double> tice = { -1., -1, -1 }; //˚C
    double hice = 0.1; // m
    double cice = 0.5;
    double dml = 10.; // m

    ModuleLoader::getLoader().setAllDefaults();
    ConfiguredModule::parseConfigurator();

    ElementData data(3);
    data.configure(); // Configure with the default linear freezing point

    data = PrognosticData::generate(hice, cice, sst, sss, 0., tice);
    data.setTimestep(86400.); // s. Very long TS to get below freezing

    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 0;
    data.incomingShortwave() = 0;

    NextsimPhysics nsphys;

    // Below the lower limit
    data.windSpeed() = 1.5;
    nsphys.updateDerivedData(data, data, data);
    nsphys.calculate(data, data, data);
    REQUIRE(0.00126936 == Approx(data.dragPressure()).epsilon(1e-4));

    // Wind speed dependent
    data.windSpeed() = 8;
    nsphys.updateDerivedData(data, data, data);
    nsphys.calculate(data, data, data);
    REQUIRE(0.00141407 == Approx(data.dragPressure()).epsilon(1e-4));

    // Above the upper limit
    data.windSpeed() = 23;
    nsphys.updateDerivedData(data, data, data);
    nsphys.calculate(data, data, data);
    REQUIRE(0.00253872 == Approx(data.dragPressure()).epsilon(1e-4));

}

TEST_CASE("Melting conditions", "[NextsimPhysics]")
{
    Configurator::clear();
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
    data.configure(); // Configure with the UNESCO freezing point

    data = PrognosticData::generate(hice, cice, sst, sss, hsnow, tice);
    data.setTimestep(600.); // s. Very long TS to get below freezing

    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 330;
    data.incomingShortwave() = 50;
    data.snowfall() = 0;

    data.windSpeed() = 5;

    NextsimPhysics nsphys;
    nsphys.configure();

    nsphys.updateDerivedData(data, data, data);
    nsphys.calculate(data, data, data);

    // Externally visible values (in PhysicsData)
    REQUIRE(0.12846 == Approx(data.updatedIceTrueThickness()).epsilon(1e-4));
    REQUIRE(0.01957732 == Approx(data.updatedSnowTrueThickness()).epsilon(1e-4));
    REQUIRE(0.368269 == Approx(data.updatedIceConcentration()).epsilon(1e-4));
    REQUIRE(0.0 == Approx(data.updatedIceSurfaceTemperature()).epsilon(1e-4));

    // Values used by Nextsim Physics modules
    REQUIRE(0.0 == nsphys.newIce());
    REQUIRE(-84.6156 == Approx(nsphys.QIceAtmosphere()).epsilon(1e-2));
    REQUIRE(53717.8 == Approx(nsphys.QIceOceanHeat()).epsilon(1e-2));
    REQUIRE(-7.3858e-06 == Approx(nsphys.sublimationRate()).epsilon(1e-4));
    REQUIRE(19.7013 == Approx(nsphys.QDerivativeWRTTemperature()).epsilon(1e-2));
    REQUIRE(0.0 == Approx(nsphys.totalIceFromSnow()).epsilon(1e-2));


}

TEST_CASE("Freezing conditions", "[NextsimPhysics]")
{
    Configurator::clear();
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

    double tair = -12; //˚C
    double tdew = -12; //˚C
    double pair = 100000; // Pa, slightly low pressure
    double sst = -1.75; //˚C
    double sss = 32; // PSU
    std::vector<double> tice = { -9., -9. }; //˚C
    double hice = 0.1; // m
    double cice = 0.5;
    double hsnow = 0.01; // m
    double dml = 10.; // m

    ModuleLoader::getLoader().setAllDefaults();
    ConfiguredModule::parseConfigurator();
    tryConfigure(ModuleLoader::getLoader().getImplementation<IIceAlbedo>());

    ElementData data(2);
    data.configure(); // Configure with the UNESCO freezing point

    data = PrognosticData::generate(hice, cice, sst, sss, hsnow, tice);
    data.setTimestep(600.); // s. Very long TS to get below freezing

    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 265;
    data.incomingShortwave() = 0;
    data.snowfall() = 1e-3;

    data.windSpeed() = 5;

    NextsimPhysics nsphys;
    nsphys.configure();

    nsphys.updateDerivedData(data, data, data);
    nsphys.calculate(data, data, data);

    // Externally visible values (in PhysicsData)
    REQUIRE(0.199998 == Approx(data.updatedIceTrueThickness()).epsilon(1e-4));
    REQUIRE(0.02179357 == Approx(data.updatedSnowTrueThickness()).epsilon(1e-4));
    REQUIRE(0.5002 == Approx(data.updatedIceConcentration()).epsilon(1e-4));
    REQUIRE(-8.90443 == Approx(data.updatedIceSurfaceTemperature()).epsilon(1e-4));

    // Values used by Nextsim Physics modules
    REQUIRE(6.79707e-5 == Approx(nsphys.newIce()).epsilon(1e-2));
    REQUIRE(42.2955 == Approx(nsphys.QIceAtmosphere()).epsilon(1e-2));
    REQUIRE(73.9465 == Approx(nsphys.QIceOceanHeat()).epsilon(1e-2));
    REQUIRE(2.15132e-06 == Approx(nsphys.sublimationRate()).epsilon(1e-4));
    REQUIRE(16.7615 == Approx(nsphys.QDerivativeWRTTemperature()).epsilon(1e-2));
    REQUIRE(0.0 == Approx(nsphys.totalIceFromSnow()).epsilon(1e-2));


}
} /* namespace Nextsim */
