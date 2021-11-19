/*!
 * @file NextsimPhysics_test.cpp
 *
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/ElementData.hpp"
#include "include/NextsimPhysics.hpp"
#include "include/constants.hpp"
namespace Nextsim {

TEST_CASE("Outgoing LW (OW)", "[NextsimPhysics]")
{
    ElementData<NextsimPhysics> data;

    const double t = 280.; // kelvin
    data = PrognosticData::generate(0., 0., celsius(t), 0., 0, { 0., 0, 0 });

    // Configure as NextsimPhysics, as the only subclass of Configured.
    data.configure();

    NextsimPhysics::heatFluxOpenWater(data, data, data, data);

    double target = PhysicalConstants::sigma * t * t * t * t;

    REQUIRE(data.QLongwaveOpenWater() == target);
}

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
} /* namespace Nextsim */
