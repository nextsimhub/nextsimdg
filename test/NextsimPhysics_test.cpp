/*!
 * @file NextsimPhysics_test.cpp
 *
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "ConfiguredModule.hpp"
#include "ElementData.hpp"
#include "NextsimPhysics.hpp"
#include "constants.hpp"
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
    ElementData<NextsimPhysics> data;

    double tair = -3;
    double tdew = 0.1;
    double pair = 100000; // Slightly low pressure
    double sst = -1;
    double sss = 32; // PSU
    std::array<double, N_ICE_TEMPERATURES> tice = { -2., -2, -2 };
    double hice = 0.1;
    double cice = 0.5;

    data = PrognosticData::generate(hice, cice, sst, sss, 0., tice);
    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;

    data.updateDerivedData(data, data, data);

    REQUIRE(1.29253 == Approx(data.airDensity()).epsilon(1e-4));
    REQUIRE(0.00385326 == Approx(data.specificHumidityAir()).epsilon(1e-4));
    REQUIRE(0.00349446 == Approx(data.specificHumidityWater()).epsilon(1e-4));
    REQUIRE(0.00323958 == Approx(data.specificHumidityIce()).epsilon(1e-4));
    REQUIRE(1020.773 == Approx(data.heatCapacityWetAir()).epsilon(1e-4));
}

TEST_CASE("New ice formation", "[NextsimPhysics]")
{
    ElementData<NextsimPhysics> data;
    ExternalData exter;
    PhysicsData phys;

    std::stringstream config;
    config << "Nextsim::IFreezingPoint = Nextsim::UnescoFreezing" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    double tair = -3; //˚C
    double tdew = 0.1; //˚C
    double pair = 100000; // Pa, slightly low pressure
    double sst = -1.5; //˚C
    double sss = 32; // PSU
    std::array<double, N_ICE_TEMPERATURES> tice = { -2., -2, -2 }; //˚C
    double hice = 0.1; // m
    double cice = 0.5;
    double dml = 10.; // m

    ConfiguredModule::parseConfigurator();
    data.configure(); // Configure with the default linear freezing point

    data = PrognosticData::generate(hice, cice, sst, sss, 0., tice);
    data.setTimestep(86400.); // s. Very long TS to get below freezing

    //data = exter;
    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 0;
    data.incomingShortwave() = 0;

    //data = phys;

    data.updateDerivedData(data, data, data);
    data.calculate(data, data, data);

    // Correct for the non-NIST value of the Stefan-Boltzman constant used by old NeXtSIM. This
    // propagates linearly through to the new sea ice value in this case.
    double sbCorr = PhysicalConstants::sigma / 5.67e-8;
    REQUIRE(0.0258236 * sbCorr == Approx(data.newIce()).epsilon(1e-4));
}

TEST_CASE("Drag pressure", "[NextsimPhysics]")
{
    ElementData<NextsimPhysics> data;
    ExternalData exter;
    PhysicsData phys;

    double tair = 2; //˚C
    double tdew = 1.5; //˚C
    double pair = 100000; // Pa, slightly low pressure
    double sst = -1.5; //˚C
    double sss = 32; // PSU
    std::array<double, N_ICE_TEMPERATURES> tice = { -1., -1, -1 }; //˚C
    double hice = 0.1; // m
    double cice = 0.5;
    double dml = 10.; // m

    ConfiguredModule::parseConfigurator();
    data.configure(); // Configure with the default linear freezing point

    data = PrognosticData::generate(hice, cice, sst, sss, 0., tice);
    data.setTimestep(86400.); // s. Very long TS to get below freezing

    data = exter;
    data.airTemperature() = tair;
    data.dewPoint2m() = tdew;
    data.airPressure() = pair;
    data.mixedLayerDepth() = dml;
    data.incomingLongwave() = 0;
    data.incomingShortwave() = 0;

    data = phys;

    // Below the lower limit
    data.windSpeed() = 1.5;
    data.updateDerivedData(data, data, data);
    data.calculate(data, data, data);
    REQUIRE(0.00126936 == Approx(data.dragPressure()).epsilon(1e-4));

    // Wind speed dependent
    data.windSpeed() = 8;
    data.updateDerivedData(data, data, data);
    data.calculate(data, data, data);
    REQUIRE(0.00141407 == Approx(data.dragPressure()).epsilon(1e-4));

    // Above the upper limit
    data.windSpeed() = 23;
    data.updateDerivedData(data, data, data);
    data.calculate(data, data, data);
    REQUIRE(0.00253872 == Approx(data.dragPressure()).epsilon(1e-4));

}

} /* namespace Nextsim */
