#include "dynamics.hpp"
#include "dgvisu.hpp"
#include "stopwatch.hpp"

namespace Nextsim {
extern Nextsim::Timer GlobalTimer;

/*!
 * Sets important parameters, initializes these and that
 */
void Dynamics::BasicInit()
{
    //! Degree of the time stepping in advection
    dgtransport.settimesteppingscheme("rk3");

    //! set meshes
    dgtransport.setmesh(mesh);

    //! Init Vectors
    vx.resize_by_mesh(mesh);
    vy.resize_by_mesh(mesh);
    A.resize_by_mesh(mesh);
    H.resize_by_mesh(mesh);
    S11.resize_by_mesh(mesh);
    S12.resize_by_mesh(mesh);
    S22.resize_by_mesh(mesh);

    E11.resize_by_mesh(mesh);
    E12.resize_by_mesh(mesh);
    E22.resize_by_mesh(mesh);
    D.resize_by_mesh(mesh);

    oceanX.resize_by_mesh(mesh);
    oceanY.resize_by_mesh(mesh);
    atmX.resize_by_mesh(mesh);
    atmY.resize_by_mesh(mesh);

    tmpX.resize_by_mesh(mesh);
    tmpY.resize_by_mesh(mesh);
}

//////////////////////////////////////////////////

/**!
 * controls the flow of the dynamical core
 */

void Dynamics::advectionStep(const double dt)
{
    GlobalTimer.start("dyn -- adv -- reinit");
    dgtransport.reinitvelocity();
    GlobalTimer.stop("dyn -- adv -- reinit");

    GlobalTimer.start("dyn -- adv -- step");
    dgtransport.step(dt, A); // performs one time step with the 2nd Order Heun scheme
    dgtransport.step(dt, H); // performs one time step with the 2nd Order Heun scheme

    // Limit H and A
#pragma omp parallel for
    for (size_t i = 0; i < mesh.n; ++i) {
        A(i, 0) = std::max(std::min(A(i, 0), 1.0), 0.0);
        H(i, 0) = std::max(H(i, 0), 0.0);
    }

    GlobalTimer.stop("dyn -- adv -- step");
}

void Dynamics::velocityContinuity(double gamma)
{
    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            velocityContinuityY(ic, ic + 1, gamma);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            velocityContinuityX(ic, ic + mesh.nx, gamma);
    }
}

void Dynamics::addStressTensorEdges(double scaleSigma)
{

    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            stressTensorEdgesY(scaleSigma, ic, ic + 1);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            stressTensorEdgesX(scaleSigma, ic, ic + mesh.nx);
    }
}

void Dynamics::addStressTensorBoundary(double scaleSigma)
{
    //! consistency term on lower and upper boundary
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {

        const size_t clower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;

        stressTensorBoundaryUpper(scaleSigma, cupper);
        stressTensorBoundaryLower(scaleSigma, clower);
    }

    //! consistency term on left and right boundary
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t cright = (iy + 1) * mesh.nx - 1;

        stressTensorBoundaryLeft(scaleSigma, cleft);
        stressTensorBoundaryRight(scaleSigma, cright);
    }
}

void Dynamics::momentumSymmetry() { }

void Dynamics::velocityDirichletBoundary(double gamma)
{

    for (size_t ix = 0; ix < mesh.nx; ++ix) {

        const size_t clower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;

        velocityDirichletBoundaryTop(cupper, gamma);
        velocityDirichletBoundaryBottom(clower, gamma);
    }

    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t cright = (iy + 1) * mesh.nx - 1;

        velocityDirichletBoundaryLeft(cleft, gamma);
        velocityDirichletBoundaryRight(cright, gamma);
    }
}

void Dynamics::computeStrainRateTensor()
{
    E11.col(0) = 1. / mesh.hx * vx.col(1);
    E11.col(1) = 1. / mesh.hx * 2. * vx.col(3);
    E11.col(2) = 1. / mesh.hx * vx.col(5);
    E12.col(0) = 1. * 0.5 * (vy.col(1) / mesh.hx + vx.col(2) / mesh.hy);
    E12.col(1) = 1. * 0.5 * (vx.col(5) / mesh.hy + 2.0 * vy.col(3) / mesh.hx);
    E12.col(2) = 1. * 0.5 * (2.0 * vx.col(4) / mesh.hy + vy.col(5) / mesh.hx);
    E22.col(0) = 1. / mesh.hy * vy.col(2);
    E22.col(1) = 1. / mesh.hy * vy.col(5);
    E22.col(2) = 1. / mesh.hy * 2. * vy.col(4);
}

//! Computes the cell terms of sigma = 1/2(nabla v + nabla v^T)
void Dynamics::addStressTensorCell(double scaleSigma)
{
    // S11 d_x phi_x + S12 d_y phi_x
    tmpX.col(1) += 12. * scaleSigma / mesh.hx * (S11.col(0));
    tmpX.col(3) += 180. * scaleSigma / mesh.hx * (S11.col(1) / 6.);
    tmpX.col(5) += 144. * scaleSigma / mesh.hx * (S11.col(2) / 12.);

    tmpX.col(2) += 12. * scaleSigma / mesh.hy * (S12.col(0));
    tmpX.col(4) += 180. * scaleSigma / mesh.hy * (S12.col(2) / 6.);
    tmpX.col(5) += 144. * scaleSigma / mesh.hy * (S12.col(1) / 12.);

    // S12 d_x phi_y + S22 d_y phi_y
    tmpY.col(1) += 12. * scaleSigma / mesh.hx * (S12.col(0));
    tmpY.col(3) += 180. * scaleSigma / mesh.hx * (S12.col(1) / 6.);
    tmpY.col(5) += 144. * scaleSigma / mesh.hx * (S12.col(2) / 12.);

    tmpY.col(2) += 12. * scaleSigma / mesh.hy * (S22.col(0));
    tmpY.col(4) += 180. * scaleSigma / mesh.hy * (S22.col(2) / 6.);
    tmpY.col(5) += 144. * scaleSigma / mesh.hy * (S22.col(1) / 12.);
}

void Dynamics::addStressTensor(double scaleSigma)
{
    addStressTensorCell(scaleSigma); //!< cell term
    addStressTensorEdges(scaleSigma);
    addStressTensorBoundary(scaleSigma);
}

}