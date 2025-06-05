/*
 * @file SpecificHumidity_test.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/FiniteElementSpecHum.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("FiniteElementSpecHum");
TEST_CASE("Specific humidity test")
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
    REQUIRE(0.00385326 == doctest::Approx(air).epsilon(prec));
    REQUIRE(0.00349446 == doctest::Approx(water).epsilon(prec));
    REQUIRE(0.00323958 == doctest::Approx(ice).epsilon(prec));
}
TEST_SUITE_END();

} /* namespace Nextsim */
