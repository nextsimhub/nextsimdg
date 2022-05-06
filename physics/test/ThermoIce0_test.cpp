/*!
 * @file ThermoIce0_test.cpp
 *
 * @date Sep 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoIce0ED.hpp"

#include <memory>

#include "include/PrognosticElementData.hpp"
#include "include/ExternalData.hpp"
#include "include/NextsimPhysics.hpp"
#include "include/PhysicsData.hpp"
#include "include/constants.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

namespace Nextsim {

TEST_CASE("Test no ice", "[ThermoIce0]")
{
    PrognosticElementData prog;
    ExternalData exter;
    PhysicsData phys;
    NextsimPhysics nsp;

    ThermoIce0ED ti0;
    double freezingPointIce = -Water::mu * Ice::s;

    prog = PrognosticElementData::generate(0, 0.99, 25, 34.56, 0, { -0.5, -0.6, -0.7 });
    phys.iceTrueThickness() = 0.25; // An arbitrary non-zero value
    phys.snowTrueThickness() = 0.10; // An arbitrary non-zero value
    phys.updatedIceSurfaceTemperature() = 0;

    ti0.calculate(prog, exter, phys, nsp);

    REQUIRE(phys.iceTrueThickness() == 0);
    REQUIRE(phys.snowTrueThickness() == 0);
    REQUIRE(phys.updatedIceSurfaceTemperature() == freezingPointIce);

    prog = PrognosticElementData::generate(0.25, 0, 25, 34.56, 0, { -0.5, -0.6, -0.7 });
    phys.iceTrueThickness() = 0.25; // An arbitrary non-zero value
    phys.snowTrueThickness() = 0.10; // An arbitrary non-zero value
    phys.updatedIceSurfaceTemperature() = 0;

    ti0.calculate(prog, exter, phys, nsp);

    REQUIRE(phys.iceTrueThickness() == 0);
    REQUIRE(phys.snowTrueThickness() == 0);
    REQUIRE(phys.updatedIceSurfaceTemperature() == freezingPointIce);
}
} /* namespace Nextsim */
