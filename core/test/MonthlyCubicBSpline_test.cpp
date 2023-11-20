/*!
 * @file MonthlyCubicBSpline.hpp
 *
 * @date Nov  7, 2023
 * @author Einar Örn Ólason <einar.olason@nersc.no>
 */

#include "include/MonthlyCubicBSpline.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

namespace Nextsim {

/* Tests for the implementation of a cyclic, monthly cubic b-spline. It's based on boost's b-spline,
 * so we only need to test if this implementation is truly cyclic and if the mapping to months works
 * as expected. The b-spline itself is assumed to work. */
TEST_SUITE_BEGIN("MonthlyCubicBSpline");

/* Test whether we get a cyclic results. This is done by giving the function a non-trivial set of
 * inputs and test if we get the same results for day 0 as we get for day 365 (or 366, if it's a
 * leap year). */
TEST_CASE("Test cyclicity")
{
    const std::vector<double> testTable = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12. };

    monthlyCubicBSpline splineTest(testTable);

    const double prec = 1e-6;
    REQUIRE(splineTest(0, 0) == doctest::Approx(splineTest(365., 0)).epsilon(prec));
    REQUIRE(splineTest(0, 1) == doctest::Approx(splineTest(366., 1)).epsilon(prec));
}

/* Test whether we get approximately right values at mid-month for a realistic set of input values.
 * Here we use the albedo values from Maykut and Untersteiner (1971). We will never get exactly the
 * tabulated value from a b-spline (except for some special cases), so we need to decide what is
 * acceptable accuracy. It is set to prec = 5e-4 here for tabulated values accurate to O(1e-2). */
TEST_CASE("Albedo values")
{
    // Monthly snow albedo from Maykut and Untersteiner (1971)
    const std::vector<double> albedoTable
        = { 0.85, 0.85, 0.83, 0.81, 0.82, 0.78, 0.64, 0.69, 0.84, 0.85, 0.85, 0.85 };
    //      Jan,  Feb,  Mar,  Apr,  Mai,  Jun,  Jul,  Aug,  Sept, Oct,  Nov,  Dec

    monthlyCubicBSpline splineTest(albedoTable);

    const double prec = 5e-4;

    // Loop over the year by mid-month days-of-year (which can be a fraction of a day).
    // Start on January 14th at 5:14:33 and assume each month is 1/12th of a Gregorian year.
    double dayOfYear = 365.2425 / 24;
    for (auto& element : albedoTable) {
        REQUIRE(splineTest(dayOfYear, 0) == doctest::Approx(element).epsilon(prec));
        dayOfYear += 365.2425 / 12;
    }
}
}
