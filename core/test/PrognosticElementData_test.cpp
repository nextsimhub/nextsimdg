/*!
 * @file PrognosticElementData_test.cpp
 * @date Dec 22, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/include/PrognosticElementData.hpp"
#include "include/PrognosticGenerator.hpp"

namespace Nextsim {

TEST_CASE("Ice-layer access function", "[PrognosticElementData]")
{
    std::vector<double> tice = { -0.1, -0.2, -0.3 };
    PrognosticElementData pd(
        PrognosticGenerator().hice(0.1).cice(0.5).hsnow(0.).tice(tice).sst(-1.).sss(32.));
    tryConfigure(pd);

    REQUIRE(pd.iceTemperature(0) == tice[0]);
    REQUIRE(pd.iceTemperature(1) == tice[1]);
    REQUIRE(pd.iceTemperature(2) == tice[2]);
}

} /* namespace Nextsim */
