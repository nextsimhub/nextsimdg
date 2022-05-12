/*!
 * @file StructureFactory_test.cpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/StructureFactory.hpp"
#include "include/DevGrid.hpp"
#include "include/DevGridIO.hpp"

#include <cstdio>

namespace Nextsim {

TEST_CASE("A valid structure name", "[StructureFactory]")
{
    DevGrid grid;
}

TEST_CASE("An invalid structure name", "[StructureFactory]")
{
}

TEST_CASE("Read a structure name from file", "[StructureFactory]")
{
}

}
