/*!
 * @file PrognosticData_test.cpp
 * @date Dec 22, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/PrognosticData.hpp"
#include "include/ModuleLoader.hpp"

namespace Nextsim {

TEST_CASE("Ice-layer access function", "[PrognosticData]")
{
    std::vector<double> tice = { -0.1, -0.2, -0.3 };
    PrognosticData pd(3);
    pd = PrognosticData::generate(0.1, 0.5, -1., 32, 0, tice);
    ModuleLoader::getLoader().setAllDefaults();
    tryConfigure(pd);

    REQUIRE(pd.iceTemperature(0) == tice[0]);
    REQUIRE(pd.iceTemperature(1) == tice[1]);
    REQUIRE(pd.iceTemperature(2) == tice[2]);
}

} /* namespace Nextsim */
