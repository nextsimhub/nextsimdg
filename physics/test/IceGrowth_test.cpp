/*!
 * @file IceGrowth_test.cpp
 *
 * @date Apr 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <sstream>

#include "include/ConfiguredModule.hpp"
#include "include/ModelArray.hpp"

namespace Nextsim {

TEST_CASE("New ice formation", "[IceGrowth]")
{
    ModelArray::setDimensions(ModelArray::Type::H, {1});
}

}
