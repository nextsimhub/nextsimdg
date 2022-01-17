#include "newdynamics.hpp"
#include "dgvisu.hpp"
#include "stopwatch.hpp"


namespace Nextsim {
extern Nextsim::Timer GlobalTimer;

void NewDynamics::velocityContinuity(double gamma, const Mesh& mesh, CellVector<2>& tmpX,
    CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
{
    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            velocityContinuityY(ic, ic + 1, gamma, mesh, tmpX, tmpY, vx, vy);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            velocityContinuityX(ic, ic + mesh.nx, gamma, mesh, tmpX, tmpY, vx, vy);
    }
}

void NewDynamics::addStressTensorEdges(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
    CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
{

    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            stressTensorEdgesY(scaleSigma, ic, ic + 1, mesh, tmpX, tmpY, S11, S12, S22);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            stressTensorEdgesX(scaleSigma, ic, ic + mesh.nx, mesh, tmpX, tmpY, S11, S12, S22);
    }
}

void NewDynamics::addStressTensorBoundary(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
    CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
{
    //! consistency term on lower and upper boundary
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {

        const size_t clower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;

        stressTensorBoundaryUpper(scaleSigma, cupper, mesh, tmpX, tmpY, S11, S12, S22);
        stressTensorBoundaryLower(scaleSigma, clower, mesh, tmpX, tmpY, S11, S12, S22);
    }

    //! consistency term on left and right boundary
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t cright = (iy + 1) * mesh.nx - 1;

        stressTensorBoundaryLeft(scaleSigma, cleft, mesh, tmpX, tmpY, S11, S12, S22);
        stressTensorBoundaryRight(scaleSigma, cright, mesh, tmpX, tmpY, S11, S12, S22);
    }
}

void NewDynamics::momentumSymmetry() { }

void NewDynamics::velocityDirichletBoundary(double gamma, const Mesh& mesh, CellVector<2>& tmpX,
    CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
{

    for (size_t ix = 0; ix < mesh.nx; ++ix) {

        const size_t clower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;

        velocityDirichletBoundaryTop(cupper, gamma, mesh, tmpX, tmpY, vx, vy);
        velocityDirichletBoundaryBottom(clower, gamma, mesh, tmpX, tmpY, vx, vy);
    }

    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t cright = (iy + 1) * mesh.nx - 1;

        velocityDirichletBoundaryLeft(cleft, gamma, mesh, tmpX, tmpY, vx, vy);
        velocityDirichletBoundaryRight(cright, gamma, mesh, tmpX, tmpY, vx, vy);
    }
}

void NewDynamics::computeStrainRateTensor(const Mesh& mesh, CellVector<2>& vx,
    CellVector<2>& vy, CellVector<1>& E11, CellVector<1>& E12, CellVector<1>& E22)
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
void NewDynamics::addStressTensorCell(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
    CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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

void NewDynamics::addStressTensor(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
    CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
{
    addStressTensorCell(scaleSigma, mesh, tmpX, tmpY, S11, S12, S22); //!< cell term
    addStressTensorEdges(scaleSigma, mesh, tmpX, tmpY, S11, S12, S22);
    addStressTensorBoundary(scaleSigma, mesh, tmpX, tmpY, S11, S12, S22);
}

}
