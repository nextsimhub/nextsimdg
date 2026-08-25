/*
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/FiniteElementSpecHum.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("FiniteElementSpecHum");
TEST_CASE("Specific humidity test")
{
    //   FloatType tair = -3;
    FloatType tdew = 0.1;
    FloatType pair = 100000; // Slightly low pressure
    FloatType sst = -1;
    FloatType sss = 32; // PSU
    std::vector<FloatType> tice = { -2., -2, -2 };

    FiniteElementSpecHum& feshw = FiniteElementSpecHum::water();
    FiniteElementSpecHum& feshi = FiniteElementSpecHum::ice();
    FloatType water = feshw(sst, pair, sss);
    FloatType air = feshw(tdew, pair);
    FloatType ice = feshi(tice[0], pair);

    FloatType prec = 1e-4;
    REQUIRE(0.00385326 == doctest::Approx(air).epsilon(prec));
    REQUIRE(0.00349446 == doctest::Approx(water).epsilon(prec));
    REQUIRE(0.00323958 == doctest::Approx(ice).epsilon(prec));
}
TEST_SUITE_END();

} /* namespace Nextsim */
