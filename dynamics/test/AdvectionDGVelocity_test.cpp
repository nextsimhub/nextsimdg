/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @brief Test the transport with a DG velocity (DGTransport::prepareAdvectionDG), as used by
 * the LDG momentum.
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "DGTransport.hpp"
#include "Interpolations.hpp"
#include "ParametricMesh.hpp"

#include <cmath>

namespace Nextsim {

static const FloatType L = 1000.0;

class Bump : public Interpolations::Function {
public:
    FloatType operator()(FloatType x, FloatType y) const
    {
        const FloatType r2 = (std::pow(x / L - 0.4, 2) + std::pow(y / L - 0.5, 2)) / 0.04;
        return (r2 < 1.0) ? std::exp(-1.0 / (1.0 - r2)) : 0.0;
    }
};
class RotationX : public Interpolations::Function {
public:
    FloatType operator()(FloatType x, FloatType y) const { return -(y - 0.5 * L) / L; }
};
class RotationY : public Interpolations::Function {
public:
    FloatType operator()(FloatType x, FloatType y) const { return (x - 0.5 * L) / L; }
};

//! Cartesian mesh of [0,L]^2 with Nx x Ny cells, optionally with an outer layer of land
static void createMesh(ParametricMesh& smesh, size_t Nx, size_t Ny, bool landFrame = true)
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
            smesh.vertices(ii, 0) = L * ix / Nx;
            smesh.vertices(ii, 1) = L * iy / Ny;
        }
    smesh.landmask.assign(smesh.nelements, 1);
    for (size_t iy = 0; iy < Ny; ++iy)
        for (size_t ix = 0; ix < Nx; ++ix)
            if (landFrame && (ix == 0 || iy == 0 || ix + 1 == Nx || iy + 1 == Ny))
                smesh.landmask[iy * Nx + ix] = 0;
}

//! A constant velocity is integrated exactly by both paths: projected and full velocity
//! must give the same result.
template <int DG> static void compareConstantVelocity(const ParametricMesh& smesh, int nsteps = 20)
{
    constexpr int DGm = (DG < 3) ? 3 : DG;
    DGVector<DG> vx(smesh), vy(smesh);
    vx.setZero();
    vy.setZero();
    vx.col(0).setConstant(0.6);
    vy.col(0).setConstant(-0.3);
    DGVector<DGm> vxm(smesh), vym(smesh);
    vxm.setZero();
    vym.setZero();
    vxm.col(0).setConstant(0.6);
    vym.col(0).setConstant(-0.3);

    DGVector<DG> phi1(smesh), phi2(smesh);
    Interpolations::Function2DG(smesh, phi1, Bump());
    phi2 = phi1;

    DGTransport<DG> projected(smesh), full(smesh);
    projected.setVelocities(vx, vy);
    projected.reinitnormalvelocity();
    full.prepareAdvectionDG(vxm, vym);

    const FloatType mass0 = phi1.col(0).sum();
    for (int it = 0; it < nsteps; ++it) {
        projected.step(10.0, phi1);
        full.step(10.0, phi2);
    }
    CHECK((phi1 - phi2).norm() / phi1.norm() < 1.e-10);
    if (nsteps > 50) // the bump has partially left the domain
        CHECK(phi1.col(0).sum() < 0.99 * mass0);
}

//! Rotating DQ1 velocity: mass is conserved (land frame, no outflow)
template <int DG> static void checkMassRotation(const ParametricMesh& smesh)
{
    DGVector<4> vx(smesh), vy(smesh);
    Interpolations::Function2DG(smesh, vx, RotationX());
    Interpolations::Function2DG(smesh, vy, RotationY());

    DGVector<DG> phi(smesh);
    Interpolations::Function2DG(smesh, phi, Bump());
    const FloatType mass0 = phi.col(0).sum();

    DGTransport<DG> transport(smesh);
    transport.prepareAdvectionDG(vx, vy);
    for (int it = 0; it < 50; ++it)
        transport.step(5.0, phi);

    CHECK(std::fabs(phi.col(0).sum() - mass0) / mass0 < 1.e-12);
    CHECK(phi.norm() > 0.5 * mass0 / std::sqrt(smesh.nelements)); // something is left
}

TEST_SUITE_BEGIN("AdvectionDGVelocity");
TEST_CASE("Constant velocity: full and projected velocity agree")
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 20, 16);
    compareConstantVelocity<1>(smesh);
    compareConstantVelocity<3>(smesh);
    compareConstantVelocity<4>(smesh);
    compareConstantVelocity<6>(smesh);
}
TEST_CASE("Constant velocity, free outer boundary: full and projected velocity agree")
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 20, 16, false);
    compareConstantVelocity<1>(smesh, 80);
    compareConstantVelocity<3>(smesh, 80);
    compareConstantVelocity<4>(smesh, 80);
    compareConstantVelocity<6>(smesh, 80);
}
TEST_CASE("Rotating DG velocity: mass conservation")
{
    ParametricMesh smesh(CARTESIAN);
    createMesh(smesh, 20, 16);
    checkMassRotation<1>(smesh);
    checkMassRotation<3>(smesh);
    checkMassRotation<4>(smesh);
    checkMassRotation<6>(smesh);
}
TEST_SUITE_END();

}
