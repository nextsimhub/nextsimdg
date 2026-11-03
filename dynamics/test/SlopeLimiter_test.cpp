/*!
 * @brief Tests of the vertex-based slope limiter (SlopeLimiter.hpp): smooth fields are not
 * limited in the interior, the linear part of a step function is bounded by the cell means.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "SlopeLimiter.hpp"

#include <cmath>

namespace Nextsim {

static const FloatType EPS = (sizeof(FloatType) == sizeof(float)) ? 1.e-5 : 1.e-10;

//! [0,1]^2 with N x N cells
static void createMesh(ParametricMesh& smesh, size_t N)
{
    smesh.reset();
    smesh.nx = N;
    smesh.ny = N;
    smesh.nnodes = (N + 1) * (N + 1);
    smesh.nelements = N * N;
    smesh.vertices.resize(smesh.nnodes, 2);
    size_t ii = 0;
    for (size_t iy = 0; iy <= N; ++iy)
        for (size_t ix = 0; ix <= N; ++ix, ++ii) {
            smesh.vertices(ii, 0) = FloatType(ix) / N;
            smesh.vertices(ii, 1) = FloatType(iy) / N;
        }
    smesh.landmask.assign(smesh.nelements, 1);
}

struct Linear : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return 0.3 + 2.0 * x - y; }
};
//! the derivatives are bilinear (for dQ2 with the term Y(X^2-1/12))
struct Quadratic : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const
    {
        return x * x + 0.5 * x * y - y * y + 0.3 * x * x * y - 0.2 * x * y * y;
    }
};
struct Step : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return (x + 0.3 * y < 0.55) ? 1. : 0.; }
};

//! a smooth field is not changed in the cells away from the boundary
template <int DG> static void checkSmoothNotLimited(const Interpolations::Function& f)
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 12);
    DGVector<DG> phi(smesh);
    Interpolations::Function2DG(smesh, phi, f);
    const DGVector<DG> phi0 = phi;

    SlopeLimiter<DG> SL(smesh);
    SL.limit(phi);

    FloatType maxdiff = 0;
    for (size_t iy = 1; iy + 1 < smesh.ny; ++iy)
        for (size_t ix = 1; ix + 1 < smesh.nx; ++ix) {
            const size_t c = iy * smesh.nx + ix;
            maxdiff = std::max(maxdiff, (phi.row(c) - phi0.row(c)).norm() / phi0.row(c).norm());
        }
    CHECK(maxdiff < EPS);
}

//! step function: the linear (dQ1: bilinear) part of the limited field is bounded by the cell
//! means in the corners
template <int DG> static void checkStepBounded()
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 16);
    DGVector<DG> phi(smesh);
    Interpolations::Function2DG(smesh, phi, Step());
    const FloatType minMean = phi.col(0).minCoeff();
    const FloatType maxMean = phi.col(0).maxCoeff();

    SlopeLimiter<DG> SL(smesh);
    SL.limit(phi);

    bool limited = false;
    for (size_t c = 0; c < smesh.nelements; ++c)
        for (int sx = -1; sx <= 1; sx += 2)
            for (int sy = -1; sy <= 1; sy += 2) {
                FloatType v = phi(c, 0) + 0.5 * sx * phi(c, 1) + 0.5 * sy * phi(c, 2);
                // dQ1: the bilinear term is limited together with the linear part
                if constexpr (DG == 4)
                    v += 0.25 * sx * sy * phi(c, 3);
                CHECK(v >= minMean - EPS);
                CHECK(v <= maxMean + EPS);
            }
    // the projection of the step overshoots, the limiter must have acted
    DGVector<DG> phi0(smesh);
    Interpolations::Function2DG(smesh, phi0, Step());
    for (size_t c = 0; c < smesh.nelements; ++c)
        limited = limited || ((phi.row(c) - phi0.row(c)).norm() > 1.e-6);
    CHECK(limited);
}

TEST_SUITE_BEGIN("SlopeLimiter");
TEST_CASE("Linear field is not limited")
{
    checkSmoothNotLimited<3>(Linear());
    checkSmoothNotLimited<4>(Linear());
    checkSmoothNotLimited<6>(Linear());
    checkSmoothNotLimited<9>(Linear());
}
TEST_CASE("Field with bilinear derivatives is not limited")
{
    checkSmoothNotLimited<9>(Quadratic());
}
TEST_CASE("Step function is bounded")
{
    checkStepBounded<3>();
    checkStepBounded<4>();
    checkStepBounded<6>();
    checkStepBounded<9>();
}
TEST_SUITE_END();

}
