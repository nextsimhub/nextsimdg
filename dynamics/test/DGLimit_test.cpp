/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @brief Test the limiters LimitMin / LimitMax of dgLimit.hpp
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "dgVector.hpp"

#include "codeGenerationDGinGauss.hpp"
#include "dgLimit.hpp"

namespace Nextsim {

static const FloatType EPS = (sizeof(FloatType) == sizeof(float)) ? 1.e-5 : 1.e-12;

//! values of the dg function in the 3x3 Gauss points and the 4 corners
template <int DG> static Eigen::Matrix<FloatType, 1, 13> values(const DGVector<DG>& v, size_t i)
{
    Eigen::Matrix<FloatType, 1, 13> res;
    res << v.row(i) * PSI<DG, 3>, v.row(i) * PSILagrange<DG, 2>;
    return res;
}

template <int DG> static void checkLimit()
{
    const size_t n = 50;
    DGVector<DG> v;
    v.resize(n, DG);
    for (size_t i = 0; i < n; ++i)
        for (int j = 0; j < DG; ++j)
            v(i, j) = (j == 0) ? 0.5 + 0.4 * std::sin(1.7 * i) : 0.8 * std::sin(0.9 * i + 2.1 * j);

    DGVector<DG> vmin = v, vmax = v;
    LimitMin(vmin, 0.0);
    LimitMax(vmax, 1.0);

    for (size_t i = 0; i < n; ++i) {
        // the mean value is kept (it is within the bounds here)
        CHECK(vmin(i, 0) == doctest::Approx(v(i, 0)));
        CHECK(vmax(i, 0) == doctest::Approx(v(i, 0)));
        // DG3 is limited conservatively in the corners, DG4 exactly in the corners, DG6 in the
        // Gauss points
        if (DG == 6) {
            CHECK((vmin.row(i) * PSI<DG, 3>).minCoeff() >= -EPS);
            CHECK((vmax.row(i) * PSI<DG, 3>).maxCoeff() <= 1.0 + EPS);
        } else {
            CHECK(values(vmin, i).minCoeff() >= -EPS);
            CHECK(values(vmax, i).maxCoeff() <= 1.0 + EPS);
        }
        // the shape is only scaled, never flipped
        if (v.row(i).tail(DG - 1).norm() > 0) {
            CHECK(vmin.row(i).tail(DG - 1).dot(v.row(i).tail(DG - 1)) >= 0.0);
            CHECK(vmax.row(i).tail(DG - 1).dot(v.row(i).tail(DG - 1)) >= 0.0);
        }
    }
}

TEST_SUITE_BEGIN("DGLimit");
TEST_CASE("LimitMin and LimitMax")
{
    checkLimit<3>();
    checkLimit<4>();
    checkLimit<6>();
}
TEST_SUITE_END();

}
