/*----------------------------   initial.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __initial_H
#define __initial_H
/*----------------------------   initial.h     ---------------------------*/

#include "dgvector.hpp"
#include "mesh.hpp"

namespace Nextsim {

class InitialBase {
public:
    virtual double operator()(double x, double y) const = 0;
};

class SmoothInitial : virtual public InitialBase {
public:
    virtual double operator()(double x, double y) const
    {
        double r2 = ((x - 0.7) * (x - 0.7) + (y - 0.7) * (y - 0.7));
        return exp(-200.0 * r2);
    }
};

class PyramidInitial : virtual public InitialBase {
public:
    virtual double operator()(double x, double y) const
    {
        double r2 = ((x - 0.3) * (x - 0.3) + (y - 0.6) * (y - 0.6));
        return std::max(0.0, 1.0 - 10.0 * sqrt(r2));
    }
};

class BoxInitial : virtual public InitialBase {
public:
    virtual double operator()(double x, double y) const
    {
        double r2 = ((x - 0.65) * (x - 0.65) + (y - 0.3) * (y - 0.3));
        if (sqrt(r2) < 0.1)
            return 1.0;
        return 0.0;
    }
};

class MixedInitial : virtual public InitialBase {
    virtual double operator()(double x, double y) const
    {
        //    return BoxInitial()(x, y);
        return SmoothInitial()(x, y) + BoxInitial()(x, y) + PyramidInitial()(x, y);
    }
};

//////////////////////////////////////////////////

template <int DGdegree>
void L2ProjectInitial(const Mesh& mesh,
    CellVector<DGdegree>& phi,
    const InitialBase& initial);

template <>
void L2ProjectInitial(const Mesh& mesh,
    CellVector<0>& phi,
    const InitialBase& initial)
{
    phi.setZero();

    size_t ic = 0;
    for (size_t iy = 0; iy < mesh.ny; ++iy)
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            Vertex xm = mesh.vertex(ix, iy);
            phi(ic, 0) = initial(xm[0], xm[1]);
        }
}

template <>
void L2ProjectInitial(const Mesh& mesh,
    CellVector<1>& phi,
    const InitialBase& initial)
{
    phi.setZero();

    std::array<double, 3> imass = { 1.0, 12.0, 12.0 }; // without 'h'
    std::array<double, 2> g2 = { -1.0 / sqrt(12.0), 1.0 / sqrt(12.0) };

    size_t ic = 0;
    for (size_t iy = 0; iy < mesh.ny; ++iy)
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            Vertex xm = mesh.vertex(ix, iy);
            for (unsigned short int gx = 0; gx < 2; ++gx)
                for (unsigned short int gy = 0; gy < 2; ++gy) {
                    // (f, phi0)
                    phi(ic, 0) += imass[0] * 0.25 * initial(xm[0] + mesh.hx * g2[gx], xm[1] + mesh.hx * g2[gy]);
                    phi(ic, 1) += imass[1] * 0.25 * initial(xm[0] + mesh.hx * g2[gx], xm[1] + mesh.hx * g2[gy]) * g2[gx];
                    phi(ic, 2) += imass[2] * 0.25 * initial(xm[0] + mesh.hx * g2[gx], xm[1] + mesh.hx * g2[gy]) * g2[gy];
                }
        }
}

template <>
void L2ProjectInitial(const Mesh& mesh,
    CellVector<2>& phi,
    const InitialBase& initial)
{
    phi.setZero();

    std::array<double, 6> imass = { 1.0, 12.0, 12.0, 180., 180., 144. }; // without 'h'
    std::array<double, 3> g3 = { -sqrt(3.0 / 5.0) * 0.5, 0.0, sqrt(3.0 / 5.0) * 0.5 };
    std::array<double, 3> w3 = { 5. / 18., 8. / 18., 5. / 18. };

    size_t ic = 0;
    for (size_t iy = 0; iy < mesh.ny; ++iy)
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++ic) {
            Vertex xm = mesh.vertex(ix, iy);
            for (unsigned short int gx = 0; gx < 3; ++gx)
                for (unsigned short int gy = 0; gy < 3; ++gy) {
                    double x = xm[0] + mesh.hx * g3[gx];
                    double y = xm[1] + mesh.hy * g3[gy];

                    double X = 0.5 + g3[gx];
                    double Y = 0.5 + g3[gy];
                    // (f, phi0)
                    phi(ic, 0) += imass[0] * w3[gx] * w3[gy] * initial(x, y) * 1.0; // 1
                    phi(ic, 1) += imass[1] * w3[gx] * w3[gy] * initial(x, y) * (X - 0.5);
                    phi(ic, 2) += imass[2] * w3[gx] * w3[gy] * initial(x, y) * (Y - 0.5);
                    phi(ic, 3) += imass[3] * w3[gx] * w3[gy] * initial(x, y) * ((X - 0.5) * (X - 0.5) - 1. / 12.);
                    phi(ic, 4) += imass[4] * w3[gx] * w3[gy] * initial(x, y) * ((Y - 0.5) * (Y - 0.5) - 1. / 12.);
                    phi(ic, 5) += imass[5] * w3[gx] * w3[gy] * initial(x, y) * (X - 0.5) * (Y - 0.5);
                }
        }
}

} // namespace Nextsim

/*----------------------------   initial.h     ---------------------------*/
/* end of #ifndef __initial_H */
#endif
/*----------------------------   initial.h     ---------------------------*/
