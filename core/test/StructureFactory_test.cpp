/*!
 * @file StructureFactory_test.cpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/StructureFactory.hpp"

namespace Nextsim {

TEST_CASE("A valid structure name", "[StructureFactory]")
{
    DevGrid grid;
    std::shared_ptr<IStructure> ps = StructureFactory::generate(grid.structureType());
    REQUIRE(ps->structureType() == grid.structureType());
}

TEST_CASE("An invalid structure name", "[StructureFactory]")
{
    REQUIRE_THROWS_AS(StructureFactory::generate("Øresundbro"), std::invalid_argument);
}

}
