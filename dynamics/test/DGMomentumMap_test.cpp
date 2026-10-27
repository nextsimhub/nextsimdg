/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @brief Test the precomputed matrices of the LDG momentum (ParametricDGMomentumMap) and the
 * matrix-free PCG solver.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "ParametricMap.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "dgPCG.hpp"

namespace Nextsim {

static const double TOL = (sizeof(FloatType) == sizeof(float)) ? 1.e-4 : 1.e-10;

//! uniform Cartesian mesh of Nx x Ny cells of size h, no land, Dirichlet on all four sides
static void createMesh(ParametricMesh& smesh, size_t Nx, size_t Ny, FloatType h)
{
    smesh.reset();
    smesh.nx = Nx;
    smesh.ny = Ny;
    smesh.nnodes = (Nx + 1) * (Ny + 1);
    smesh.nelements = Nx * Ny;

    smesh.vertices.resize(smesh.nnodes, 2);
    size_t ii = 0;
    for (size_t iy = 0; iy <= Ny; ++iy)
        for (size_t ix = 0; ix <= Nx; ++ix, ++ii) {
            smesh.vertices(ii, 0) = h * ix;
            smesh.vertices(ii, 1) = h * iy;
        }

    smesh.landmask.assign(smesh.nelements, 1);
    smesh.dirichletFromLandmask(true);
}

//! deterministic, non-trivial coefficients
template <int DG> static DGVector<DG> testVector(const ParametricMesh& smesh, FloatType seed)
{
    DGVector<DG> v(smesh);
    for (size_t c = 0; c < smesh.nelements; ++c)
        for (int j = 0; j < DG; ++j)
            v(c, j) = std::sin(seed + 1.3 * c + 0.7 * j);
    return v;
}

template <int DG> static void checkMap()
{
    const FloatType h = 2.0;
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 6, 5, h);

    ParametricDGMomentumMap<DG, DG, DG> pmap(smesh);
    pmap.InitializeDGDivMatrices();
    pmap.InitializeIMJwPSI();
    pmap.InitializeJwPSI();
    pmap.InitializeJumpPenalty(0.0, 10.0);

    DGVector<1> one(smesh);
    one.setConstant(1.0);
    pmap.InitializeWeightedDGMassMatrix(one);

    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        // weight 1: the plain mass matrix
        CHECK((pmap.WeightedDGMassMatrix[eid] - ParametricTools::massMatrix<DG>(smesh, eid)).norm()
            == doctest::Approx(0.0).epsilon(TOL));

        // (1, d/dx X) = (1, d/dy Y) = h on a cell of size h, M^-1 (1, d/dx X) = 12/h
        CHECK(pmap.dgDivS1[eid](1, 0) == doctest::Approx(h).epsilon(TOL));
        CHECK(pmap.dgDivS2[eid](2, 0) == doctest::Approx(h).epsilon(TOL));
        CHECK(pmap.dgDivS1[eid](2, 0) == doctest::Approx(0.0).epsilon(TOL));
        CHECK(pmap.iMgradTauX[eid](1, 0) == doctest::Approx(12.0 / h).epsilon(TOL));
        CHECK(pmap.iMgradTauY[eid](2, 0) == doctest::Approx(12.0 / h).epsilon(TOL));
    }

    // the jump penalty is symmetric and positive definite (Dirichlet on all sides)
    const DGVector<DG> x = testVector<DG>(smesh, 0.1);
    const DGVector<DG> y = testVector<DG>(smesh, 2.3);
    DGVector<DG> Jx(smesh), Jy(smesh);
    pmap.ApplyMassPlusJumpPenalty(0.0, 1.0, x, Jx);
    pmap.ApplyMassPlusJumpPenalty(0.0, 1.0, y, Jy);
    CHECK(dgDot(x, Jy) == doctest::Approx(dgDot(y, Jx)).epsilon(TOL));
    CHECK(dgDot(x, Jx) > 0.0);

    // no jumps of a constant: J vanishes in the interior cells
    DGVector<DG> c(smesh), Jc(smesh);
    c.setZero();
    c.col(0).setConstant(1.0);
    pmap.ApplyMassPlusJumpPenalty(0.0, 1.0, c, Jc);
    for (size_t iy = 1; iy + 1 < smesh.ny; ++iy)
        for (size_t ix = 1; ix + 1 < smesh.nx; ++ix)
            CHECK(Jc.row(iy * smesh.nx + ix).norm() == doctest::Approx(0.0).epsilon(TOL));

    // PCG for (M + dt J) u = F with the block-Jacobi preconditioner
    const FloatType dt = 0.5;
    std::vector<Eigen::Matrix<FloatType, DG, DG>> iDiag(smesh.nelements);
    for (size_t eid = 0; eid < smesh.nelements; ++eid)
        iDiag[eid]
            = (pmap.WeightedDGMassMatrix[eid] + dt * pmap.JumpPenaltyDiagonalBlock[eid]).inverse();

    auto applyA = [&](const DGVector<DG>& in, DGVector<DG>& out) {
        pmap.ApplyMassPlusJumpPenalty(1.0, dt, in, out);
    };
    auto applyPrec = [&](const DGVector<DG>& in, DGVector<DG>& out) {
        for (size_t eid = 0; eid < smesh.nelements; ++eid)
            out.row(eid) = (iDiag[eid] * in.row(eid).transpose()).transpose();
    };

    const DGVector<DG> F = testVector<DG>(smesh, 4.2);
    DGVector<DG> u(smesh), Au(smesh);
    u.setZero();
    FloatType err;
    const FloatType tol = (sizeof(FloatType) == sizeof(float)) ? 1.e-5 : 1.e-10;
    const int iterations = PCG<DG>(applyA, applyPrec, F, u, tol, 500, err);
    CHECK(iterations < 500);
    applyA(u, Au);
    DGVector<DG> r(smesh);
    r = F - Au;
    CHECK(std::sqrt(dgDot(r, r) / dgDot(F, F)) < 10 * tol);
}

TEST_SUITE_BEGIN("DGMomentumMap");
TEST_CASE("Dirichlet edges from the landmask")
{
    // 3 x 2 cells, land in the lower left cell
    //   3 4 5
    //   L 1 2
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 3, 2, 1.0);
    smesh.landmask[0] = 0;

    smesh.dirichletFromLandmask(false);
    CHECK(smesh.dirichlet[0] == std::vector<size_t>({ 3 }));
    CHECK(smesh.dirichlet[1].empty());
    CHECK(smesh.dirichlet[2].empty());
    CHECK(smesh.dirichlet[3] == std::vector<size_t>({ 1 }));
    smesh.dirichlet[1] = { 99 }; // cleared by the next call

    smesh.dirichletFromLandmask(true);
    CHECK(smesh.dirichlet[0] == std::vector<size_t>({ 1, 2, 3 }));
    CHECK(smesh.dirichlet[1] == std::vector<size_t>({ 2, 5 }));
    CHECK(smesh.dirichlet[2] == std::vector<size_t>({ 3, 4, 5 }));
    CHECK(smesh.dirichlet[3] == std::vector<size_t>({ 1, 3 }));
}
TEST_CASE("ParametricDGMomentumMap and PCG")
{
    checkMap<3>();
    checkMap<4>();
    checkMap<6>();
    checkMap<9>();
}
TEST_SUITE_END();

}
