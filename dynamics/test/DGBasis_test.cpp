/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @brief Test the ordering of the DG basis: the generated tables and the edge traces
 * (dgEdgeTools.hpp) against an independent implementation of the basis functions.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "dgVector.hpp"

#include "codeGenerationDGinGauss.hpp"
#include "dgEdgeTools.hpp"

#include <cmath>

namespace Nextsim {

static const double TOL = (sizeof(FloatType) == sizeof(float)) ? 1.e-5 : 1.e-12;

//! The DG basis on [0,1]^2: 1, X, Y, XY, X^2-1/12, Y^2-1/12, Y(X^2-1/12), X(Y^2-1/12),
//! (X^2-1/12)(Y^2-1/12) with X = x-1/2, Y = y-1/2
static double basis(int j, double x, double y)
{
    const double X = x - 0.5;
    const double Y = y - 0.5;
    const double X2 = X * X - 1. / 12.;
    const double Y2 = Y * Y - 1. / 12.;
    const double values[9] = { 1., X, Y, X * Y, X2, Y2, Y * X2, X * Y2, X2 * Y2 };
    return values[j];
}

//! The edge basis on [0,1]: 1, T, T^2-1/12 with T = t-1/2
static double edgeBasis(int j, double t)
{
    const double T = t - 0.5;
    const double values[3] = { 1., T, T * T - 1. / 12. };
    return values[j];
}

template <int DG> static double evaluate(const DGVector<DG>& v, double x, double y)
{
    double sum = 0;
    for (int j = 0; j < DG; ++j)
        sum += v(0, j) * basis(j, x, y);
    return sum;
}

template <int ED> static double evaluateEdge(const Eigen::Matrix<FloatType, 1, ED>& e, double t)
{
    double sum = 0;
    for (int j = 0; j < ED; ++j)
        sum += e(0, j) * edgeBasis(j, t);
    return sum;
}

template <int DG> static void checkBasis()
{
    // generated values in the 3x3 Gauss points, lexicographic, x fastest
    for (int j = 0; j < DG; ++j)
        for (int gy = 0; gy < 3; ++gy)
            for (int gx = 0; gx < 3; ++gx)
                CHECK(PSI<DG, 3>(j, 3 * gy + gx)
                    == doctest::Approx(basis(j, gauss_points3[gx], gauss_points3[gy]))
                        .epsilon(TOL));

    // edge traces of a generic function of the space
    DGVector<DG> v;
    v.resize(1, DG);
    for (int j = 0; j < DG; ++j)
        v(0, j) = 1.0 + 0.37 * j - 0.05 * j * j;

    for (double t : { 0.0, 0.2, 0.5, 0.9 }) {
        CHECK(evaluateEdge(leftedgeofcell<DG>(v, 0), t)
            == doctest::Approx(evaluate(v, 0.0, t)).epsilon(TOL));
        CHECK(evaluateEdge(rightedgeofcell<DG>(v, 0), t)
            == doctest::Approx(evaluate(v, 1.0, t)).epsilon(TOL));
        CHECK(evaluateEdge(bottomedgeofcell<DG>(v, 0), t)
            == doctest::Approx(evaluate(v, t, 0.0)).epsilon(TOL));
        CHECK(evaluateEdge(topedgeofcell<DG>(v, 0), t)
            == doctest::Approx(evaluate(v, t, 1.0)).epsilon(TOL));
    }
}

TEST_SUITE_BEGIN("DGBasis");
TEST_CASE("DG basis ordering and edge traces")
{
    checkBasis<1>();
    checkBasis<3>();
    checkBasis<4>();
    checkBasis<6>();
    checkBasis<8>();
    checkBasis<9>();
}
TEST_SUITE_END();

}
