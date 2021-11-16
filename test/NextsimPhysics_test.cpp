/*!
 * @file NextsimPhysics_test.cpp
 *
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/include/ElementData.hpp"
#include "../src/include/NextsimPhysics.hpp"
#include "../src/include/constants.hpp"
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
} /* namespace Nextsim */
