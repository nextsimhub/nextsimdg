/*!
 * @file dgInitial.hpp
 * @date July 10 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGINITIAL_HPP
#define __DGINITIAL_HPP

#include "Mesh.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"

#include "ParametricTools.hpp"

namespace Nextsim {

typedef std::function<double(double, double)> InitialOp;

struct SmoothInitial {
public:
    double operator()(double x, double y) const
    {
        double r2 = ((x - 0.7) * (x - 0.7) + (y - 0.7) * (y - 0.7));
        return exp(-200.0 * r2);
    }
};

struct PyramidInitial {
public:
    double operator()(double x, double y) const
    {
        double r2 = ((x - 0.3) * (x - 0.3) + (y - 0.6) * (y - 0.6));
        return std::max(0.0, 1.0 - 10.0 * sqrt(r2));
    }
};

struct BoxInitial {
public:
    double operator()(double x, double y) const
    {
        double r2 = ((x - 0.65) * (x - 0.65) + (y - 0.3) * (y - 0.3));
        if (sqrt(r2) < 0.1)
            return 1.0;
        return 0.0;
    }
};

struct MixedInitial {
public:
    double operator()(double x, double y) const
    {
        return SmoothInitial()(x, y) + BoxInitial()(x, y) + PyramidInitial()(x, y);
    }
};

////////////////////////////////////////////////// New Interface - SasipMesh


////////////////////////////////////////////////// OLD Interfvace

//! Functions to project an analytical solution into the DG spaces
void L2ProjectInitial(const Mesh& mesh, CellVector<1>& phi, InitialOp initial)
{
    phi.setZero();

#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = mesh.nx * iy;
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            const Vertex xm = mesh.midpoint(ix, iy);
            phi(ic, 0) = initial(xm[0], xm[1]);
        }
    }
}

void L2ProjectInitial(const Mesh& mesh, CellVector<3>& phi, InitialOp initial)
{
    phi.setZero();

    const std::array<const double, 3> imass = { 1.0, 12.0, 12.0 }; // without 'h'

    // Gauss quadrature on [-1/2, 1/2]
    const std::array<const double, 2> g2 = { -1.0 / sqrt(12.0), 1.0 / sqrt(12.0) };

#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx;
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            const Vertex xm = mesh.midpoint(ix, iy);
            for (unsigned short int gx = 0; gx < 2; ++gx)
                for (unsigned short int gy = 0; gy < 2; ++gy) {
                    // (f, phi0)
                    const double tmp
                        = 0.25 * initial(xm[0] + mesh.hx * g2[gx], xm[1] + mesh.hy * g2[gy]);
                    phi(ic, 0) += imass[0] * tmp;
                    phi(ic, 1) += imass[1] * tmp * g2[gx];
                    phi(ic, 2) += imass[2] * tmp * g2[gy];
                }
        }
    }
}
void L2ProjectInitial(const Mesh& mesh, CellVector<6>& phi, InitialOp initial)
{
    phi.setZero();

    const std::array<const double, 6> imass = { 1.0, 12.0, 12.0, 180., 180., 144. }; // without 'h'

    // Gauss quadrature on [-1/2, 1/2]
    const std::array<const double, 3> g3 = { -sqrt(3.0 / 5.0) * 0.5, 0.0, sqrt(3.0 / 5.0) * 0.5 };
    const std::array<const double, 3> w3 = { 5. / 18., 8. / 18., 5. / 18. };

#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx;
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            const Vertex xm = mesh.midpoint(ix, iy);
            for (unsigned short int gx = 0; gx < 3; ++gx)
                for (unsigned short int gy = 0; gy < 3; ++gy) {
                    const double x = xm[0] + mesh.hx * g3[gx];
                    const double y = xm[1] + mesh.hy * g3[gy];

                    const double X = 0.5 + g3[gx];
                    const double Y = 0.5 + g3[gy];
                    // (f, phi0)
                    const double tmp = w3[gx] * w3[gy] * initial(x, y);
                    phi(ic, 0) += imass[0] * tmp * 1.0; // 1
                    phi(ic, 1) += imass[1] * tmp * (X - 0.5);
                    phi(ic, 2) += imass[2] * tmp * (Y - 0.5);
                    phi(ic, 3) += imass[3] * tmp * ((X - 0.5) * (X - 0.5) - 1. / 12.);
                    phi(ic, 4) += imass[4] * tmp * ((Y - 0.5) * (Y - 0.5) - 1. / 12.);
                    phi(ic, 5) += imass[5] * tmp * (X - 0.5) * (Y - 0.5);
                }
        }
    }
}

//! Functions to compute the error of a DG vector w.r.t. an analytical solution measured in L2
double L2Error(const Mesh& mesh, const CellVector<6>& phi, InitialOp ex)
{
    double res = 0.0; //!< stores the integral

    // Gauss quadrature on [-1/2, 1/2]
    const std::array<const double, 3> g3 = { -sqrt(3.0 / 5.0) * 0.5, 0.0, sqrt(3.0 / 5.0) * 0.5 };
    const std::array<const double, 3> w3 = { 5. / 18., 8. / 18., 5. / 18. };

    const double hxhy = mesh.hx * mesh.hy;

#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx;
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            const Vertex xm = mesh.midpoint(ix, iy);
            for (unsigned short int gx = 0; gx < 3; ++gx)
                for (unsigned short int gy = 0; gy < 3; ++gy) {
                    const double x = xm[0] + mesh.hx * g3[gx];
                    const double y = xm[1] + mesh.hy * g3[gy];

                    const double X = 0.5 + g3[gx];
                    const double Y = 0.5 + g3[gy];

                    const double PHI = phi(ic, 0) + phi(ic, 1) * (X - 0.5) + phi(ic, 2) * (Y - 0.5)
                        + phi(ic, 3) * ((X - 0.5) * (X - 0.5) - 1.0 / 12)
                        + phi(ic, 4) * ((Y - 0.5) * (Y - 0.5) - 1.0 / 12)
                        + phi(ic, 5) * (X - 0.5) * (Y - 0.5);

#pragma omp atomic
                    res += w3[gx] * w3[gy] * hxhy * (PHI - ex(x, y)) * (PHI - ex(x, y));
                }
        }
    }
    return sqrt(res);
}

//! Functions to project an analytical solution into the DG spaces (OLD MESH)
void InterpolateCG(const Mesh& mesh, CGVector<2>& phi, InitialOp initial)
{
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == phi.rows());

#pragma omp parallel for
    for (size_t iy = 0; iy < 2 * mesh.ny + 1; ++iy) {
        const double Y = mesh.hy * 0.5 * iy;
        size_t ii = iy * (2 * mesh.nx + 1);

        for (size_t ix = 0; ix < 2 * mesh.nx + 1; ++ix, ++ii) {
            const double X = mesh.hx * 0.5 * ix;
            phi(ii) = initial(X, Y);
        }
    }
}
void InterpolateCG(const Mesh& mesh, CGVector<1>& phi, InitialOp initial)
{
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == phi.rows());

#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny + 1; ++iy) {
        const double Y = mesh.hy * iy;
        size_t ii = iy * (mesh.nx + 1);

        for (size_t ix = 0; ix < mesh.nx + 1; ++ix, ++ii) {
            const double X = mesh.hx * ix;
            phi(ii) = initial(X, Y);
        }
    }
}

} /* namespace Nextsim */

#endif /* __DGINITIAL_H */
