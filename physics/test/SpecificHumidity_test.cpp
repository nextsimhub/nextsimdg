/*
 * @file SpecificHumidity_test.cpp
 *
 * @date May 2, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/FiniteElementSpecHum.hpp"

namespace Nextsim {

TEST_CASE("Specific humidity test", "[FiniteElementSpecHum")
{
    double tair = -3;
    double tdew = 0.1;
    double pair = 100000; // Slightly low pressure
    double sst = -1;
    double sss = 32; // PSU
    std::vector<double> tice = { -2., -2, -2 };

    FiniteElementSpecHum& feshw = FiniteElementSpecHum::water();
    FiniteElementSpecHum& feshi = FiniteElementSpecHum::ice();
    double water = feshw(sst, pair, sss);
    double air = feshw(tdew, pair);
    double ice = feshi(tice[0], pair);

    double prec = 1e-4;
    REQUIRE(0.00385326 == Approx(air).epsilon(prec));
    REQUIRE(0.00349446 == Approx(water).epsilon(prec));
    REQUIRE(0.00323958 == Approx(ice).epsilon(prec));

}
} /* namespace Nextsim */
