/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @brief Tests of the LDG momentum: strain and stress divergence operators (manufactured
 * solutions), sea surface gradient and a free drift solution of the mEVP solver.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "Interpolations.hpp"
#include "LDGMomentum.hpp"
#include "ParametricMesh.hpp"

#include <cmath>
#include <iostream>

namespace Nextsim {

static const bool isDouble = (sizeof(FloatType) == sizeof(double));

//! [0,1]^2 with N x N cells, optionally distorted and with an outer layer of land.
//! Dirichlet edges from the landmask, the outer boundary is free.
static void createMesh(ParametricMesh& smesh, size_t N, FloatType distort, bool landFrame)
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
            smesh.vertices(ii, 0)
                = FloatType(ix) / N + distort * sin(M_PI * ix / N * 3.0) * sin(M_PI * iy / N);
            smesh.vertices(ii, 1)
                = FloatType(iy) / N + distort * sin(M_PI * iy / N * 2.0) * sin(M_PI * ix / N * 2.0);
        }
    smesh.landmask.assign(smesh.nelements, 1);
    if (landFrame)
        for (size_t iy = 0; iy < N; ++iy)
            for (size_t ix = 0; ix < N; ++ix)
                if (ix == 0 || iy == 0 || ix + 1 == N || iy + 1 == N)
                    smesh.landmask[iy * N + ix] = 0;
    smesh.dirichletFromLandmask(false);
}

//! manufactured velocity (not zero on the boundary) and its strain
struct UX : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return sin(M_PI * x) * cos(M_PI * y); }
};
struct UY : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return x * x + exp(y); }
};
struct E11Exact : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const
    {
        return M_PI * cos(M_PI * x) * cos(M_PI * y);
    }
};
struct E22Exact : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return exp(y); }
};
struct E12Exact : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const
    {
        return 0.5 * (-M_PI * sin(M_PI * x) * sin(M_PI * y) + 2.0 * x);
    }
};

//! L2 error of the LDG strain of the manufactured velocity, free outer boundary
template <int DG> static FloatType strainError(size_t N, FloatType distort)
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, N, distort, false);
    LDGMomentum<DG, 1> momentum(smesh);
    Interpolations::Function2DG(smesh, momentum.vx, UX());
    Interpolations::Function2DG(smesh, momentum.vy, UY());
    momentum.projectVelocityToStrain();
    // L2ErrorFunctionDG returns the squared error
    return std::sqrt(Interpolations::L2ErrorFunctionDG(smesh, momentum.E11, E11Exact())
        + Interpolations::L2ErrorFunctionDG(smesh, momentum.E12, E12Exact())
        + Interpolations::L2ErrorFunctionDG(smesh, momentum.E22, E22Exact()));
}

template <int DG> static void checkStrainConvergence(FloatType distort, FloatType minRate)
{
    const FloatType e1 = strainError<DG>(8, distort);
    const FloatType e2 = strainError<DG>(16, distort);
    const FloatType e3 = strainError<DG>(32, distort);
    const FloatType rate = std::log2(e2 / e3);
    std::cout << "LDG strain DG" << DG << " distort " << distort << ": " << e1 << " " << e2 << " "
              << e3 << "  rate " << std::log2(e1 / e2) << " " << rate << std::endl;
    CHECK(e3 < e2);
    CHECK(rate > minRate);
}

//! divergence free linear stress
struct S11Lin : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return x; }
};
struct S12Lin : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return 2.0 * x - y; }
};
struct S22Lin : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return -2.0 * y; }
};

//! cell term + edge terms of the stress divergence vanish for a divergence free stress
//! (land frame: Dirichlet edges)
template <int DG> static void checkDivergenceFreeStress(FloatType distort)
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 12, distort, true);
    LDGMomentum<DG, 1> momentum(smesh);
    Interpolations::Function2DG(smesh, momentum.S11, S11Lin());
    Interpolations::Function2DG(smesh, momentum.S12, S12Lin());
    Interpolations::Function2DG(smesh, momentum.S22, S22Lin());

    DGVector<DG> rx(smesh), ry(smesh);
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        rx.row(eid) = (momentum.pmap.dgDivS1[eid] * momentum.S11.row(eid).transpose()
            + momentum.pmap.dgDivS2[eid] * momentum.S12.row(eid).transpose())
                          .transpose();
        ry.row(eid) = (momentum.pmap.dgDivS1[eid] * momentum.S12.row(eid).transpose()
            + momentum.pmap.dgDivS2[eid] * momentum.S22.row(eid).transpose())
                          .transpose();
        if (smesh.landmask[eid] == 0) {
            rx.row(eid).setZero();
            ry.row(eid).setZero();
        }
    }
    const FloatType reference = rx.norm() + ry.norm(); // the cell term alone
    ldgSigmaFlux<DG, DG>(smesh, 0.4, momentum.S11, momentum.S12, momentum.S22, rx, ry);
    CHECK(reference > 0.1);
    CHECK((rx.norm() + ry.norm()) / reference < (isDouble ? 1.e-12 : 1.e-5));
}

struct LinearSSH : public Interpolations::Function {
    FloatType operator()(FloatType x, FloatType y) const { return 0.3 * x - 0.7 * y; }
};

//! Gradient of a linear sea surface height (given as dG0) is exact in the interior
static void checkSeaSurfaceGradient()
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 10, 0.0, false);
    LDGMomentum<3, 1> momentum(smesh);
    DGVector<1> ssh(smesh);
    Interpolations::Function2DG(smesh, ssh, LinearSSH());
    momentum.setSeaSurfaceHeight(ssh);
    for (size_t iy = 1; iy + 1 < smesh.ny; ++iy)
        for (size_t ix = 1; ix + 1 < smesh.nx; ++ix) {
            const size_t eid = iy * smesh.nx + ix;
            for (int g = 0; g < GAUSSPOINTS(3); ++g) {
                CHECK(momentum.sshGradX(eid, g) == doctest::Approx(0.3));
                CHECK(momentum.sshGradY(eid, g) == doctest::Approx(-0.7));
            }
        }
}

/*!
 * Without ice strength and Coriolis, uniform wind and ocean at rest, free outer boundary: the
 * velocity converges to the free drift FAtm |a| a = FOcean |u| u in all cells.
 */
template <int DG> static void checkFreeDrift()
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 10, 0.05, false);
    // 100 km domain
    smesh.vertices *= 1.e5;

    LDGMomentum<DG, 3> momentum(smesh);
    VPParameters params;
    params.pStar = 0.0;
    params.fc = 0.0;
    momentum.ax.col(0).setConstant(10.0);
    momentum.ay.col(0).setConstant(-5.0);

    DGVector<3> H(smesh), A(smesh);
    H.setZero();
    A.setZero();
    H.col(0).setConstant(1.0);
    A.col(0).setConstant(1.0);

    const FloatType alpha = 300.0, beta = 300.0, dt = 600.0;
    for (int step = 0; step < 120; ++step) {
        momentum.prepareIteration(params, H, A, beta, dt);
        for (int k = 0; k < 50; ++k)
            momentum.mEVPStep(params, alpha, beta, dt);
    }
    CHECK(momentum.lastPCGIterations() < momentum.pcgMaxIterations);

    const FloatType FAtm = params.CAtm * params.rhoAtm;
    const FloatType FOcean = params.COcean * params.rhoOcean;
    const FloatType absa = std::hypot(10.0, -5.0);
    const FloatType ux = std::sqrt(FAtm / FOcean) * 10.0;
    const FloatType uy = std::sqrt(FAtm / FOcean) * -5.0;
    const FloatType tol = isDouble ? 1.e-5 : 1.e-3;
    FloatType maxdev = 0;
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
        maxdev = std::max(maxdev, std::fabs(momentum.vx(eid, 0) - ux) / std::fabs(ux));
        maxdev = std::max(maxdev, std::fabs(momentum.vy(eid, 0) - uy) / std::fabs(uy));
        maxdev = std::max(maxdev, momentum.vx.row(eid).tail(DG - 1).norm() / std::fabs(ux));
    }
    std::cout << "LDG free drift DG" << DG << ": max relative deviation " << maxdev << std::endl;
    CHECK(maxdev < tol);
}

TEST_SUITE_BEGIN("LDGMomentum");
TEST_CASE("Strain of a manufactured velocity, free boundary")
{
    checkStrainConvergence<3>(0.0, 0.8);
    checkStrainConvergence<3>(0.05, 0.8);
    checkStrainConvergence<4>(0.0, 0.8);
    checkStrainConvergence<6>(0.0, 1.8);
    checkStrainConvergence<9>(0.05, 1.8);
}
TEST_CASE("Stress divergence of a divergence free stress, Dirichlet boundary")
{
    checkDivergenceFreeStress<3>(0.0);
    checkDivergenceFreeStress<4>(0.05);
    checkDivergenceFreeStress<6>(0.05);
    checkDivergenceFreeStress<9>(0.0);
}
TEST_CASE("Gradient of the sea surface height") { checkSeaSurfaceGradient(); }
TEST_CASE("Free drift")
{
    checkFreeDrift<3>();
    checkFreeDrift<4>();
}
TEST_SUITE_END();

}
