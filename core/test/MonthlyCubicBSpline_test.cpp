//
// Created by Einar Ólason on 24/05/2023.
//

#include "include/MonthlyCubicBSpline.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

namespace Nextsim {

TEST_SUITE_BEGIN("MonthlyCubicBSpline");
TEST_CASE("Test cyclicity")
{

    const std::vector<double> testTable = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12. };

    monthlyCubicBSpline splineTest(testTable);

    const double prec = 1e-6;
    REQUIRE(splineTest(0, 0) == doctest::Approx(splineTest(365., 0)).epsilon(prec));
    REQUIRE(splineTest(0, 1) == doctest::Approx(splineTest(366., 1)).epsilon(prec));
}
}
