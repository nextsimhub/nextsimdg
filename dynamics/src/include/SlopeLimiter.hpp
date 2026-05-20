/*!
 * @file    SlopeLimiter.hpp
 * @date    Apr 2, 2024
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __SLOPELIMITER_HPP
#define __SLOPELIMITER_HPP

#include "ParametricMesh.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"

namespace Nextsim {

// Slope limiter based on D. Kuzmin: A vertex-based hierarchical slope limiter for p-adaptive
// discontinuous Galerkin methods. Journal of Computational and Applied Mathematics 2010.

template <int DG> class SlopeLimiter {
    const ParametricMesh& mesh;

    //! minimum and maximum values at the mesh nodes
    CGVector<1> minV, maxV;
    CGVector<1> dxminV, dxmaxV;
    CGVector<1> dyminV, dymaxV;

    //! alpha-values for limiting
    DGVector<1> alpha, alphaX, alphaY;

public:
    SlopeLimiter<DG>(const ParametricMesh& _mesh)
        : mesh(_mesh)
    {
        if (DG >= 3) {
            minV.resize_by_mesh(mesh);
            maxV.resize_by_mesh(mesh);
            alpha.resize_by_mesh(mesh);
        }
        if (DG >= 6) {
            dxminV.resize_by_mesh(mesh);
            dxmaxV.resize_by_mesh(mesh);
            dyminV.resize_by_mesh(mesh);
            dymaxV.resize_by_mesh(mesh);
            alphaX.resize_by_mesh(mesh);
            alphaY.resize_by_mesh(mesh);
        }
        if (DG == 8) {
            std::cerr << "No limiting for DG8" << std::endl;
            abort();
        }
    }

    // gets minimum and maximum values at all mesh vertices
    void initMinMax(
        CGVector<1>& _minv, CGVector<1>& _maxv, const DGVector<DG>& phi, size_t comp = 0)
    {
        // relative indices of the four vertices in minV/maxV
        const size_t cgindices[4] = { 0, 1, mesh.nx + 1, mesh.nx + 2 };

#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < mesh.nnodes; ++i) {
            _minv(i) = 1.e9;
            _maxv(i) = -1.e9;
        }

        // parallelization in stripes
        for (size_t p = 0; p < 2; ++p)
#pragma omp parallel for schedule(static)
            for (size_t cy = p; cy < mesh.ny; cy += 2) //!< loop over all cells of the mesh
            {
                size_t c = mesh.nx * cy; // index of first cell in row
                size_t cgi = (mesh.nx + 1) * cy; // index of vertex of cell in row
                for (size_t cx = 0; cx < mesh.nx;
                     ++cx, ++c, ++cgi) //!< loop over all cells of the mesh
                {
                    const FloatType meanvalue = phi(c, comp);
                    for (size_t j = 0; j < 4; ++j) {
                        _minv(cgi + cgindices[j]) = std::min(_minv(cgi + cgindices[j]), meanvalue);
                        _maxv(cgi + cgindices[j]) = std::max(_maxv(cgi + cgindices[j]), meanvalue);
                    }
                }
            }
    }

    // truncates the averages by min or max value
    void limitMax(DGVector<DG>& phi, FloatType max) const
    {
#pragma omp parallel for
        for (size_t c = 0; c < mesh.nelements; ++c)
            phi(c, 0) = std::min(max, phi(c, 0));
    }
    // truncates the averages by min or max value
    void limitMin(DGVector<DG>& phi, FloatType min) const
    {
#pragma omp parallel for
        for (size_t c = 0; c < mesh.nelements; ++c)
            phi(c, 0) = std::max(min, phi(c, 0));
    }

    void computeAlphas(DGVector<1>& alpha, const DGVector<DG>& phi, const CGVector<1>& _min,
        const CGVector<1>& _max)
    {
        // relative indices of the four vertices in minV/maxV
        const size_t cgindices[4] = { 0, 1, mesh.nx + 1, mesh.nx + 2 };

        assert(alpha.rows() == mesh.nelements);
#pragma omp parallel for schedule(static)
        for (size_t c = 0; c < mesh.nelements; ++c) {
            size_t cx = c % mesh.nx;
            size_t cy = c / mesh.nx;
            size_t cgi = cy * (mesh.nx + 1) + cx; // index of lower-left vertex

            // values of phi in the 4 nodes: lower-left, lower-right, upper-left, upper-right
            const Eigen::Vector<FloatType, 4> vertexvalues
                = { phi(c, 0) - 0.5 * phi(c, 1) - 0.5 * phi(c, 2),
                      phi(c, 0) + 0.5 * phi(c, 1) - 0.5 * phi(c, 2),
                      phi(c, 0) - 0.5 * phi(c, 1) + 0.5 * phi(c, 2),
                      phi(c, 0) + 0.5 * phi(c, 1) + 0.5 * phi(c, 2) };

            // value of phi in the midpoint
            const FloatType midvalue = phi(c, 0);

            FloatType al = 1.0; // the limiter
            for (size_t i = 0; i < 4; ++i) {
                FloatType dv = vertexvalues[i] - midvalue; // distance to midpoint
                if (dv > 1.e-8) {
                    assert(_max(cgi + cgindices[i]) >= midvalue);
                    al = std::min(al, std::min(1.0, (_max(cgi + cgindices[i]) - midvalue) / dv));
                }
                if (dv < -1.e-8) {
                    assert(_min(cgi + cgindices[i]) <= midvalue);
                    al = std::min(al, std::min(1.0, (_min(cgi + cgindices[i]) - midvalue) / dv));
                }
                assert(al >= 0);
            }
            alpha(c) = al;
        }
    }

    void computeAlphasX(DGVector<1>& alphax, const DGVector<DG>& phi, const CGVector<1>& _min,
        const CGVector<1>& _max)
    {
        // relative indices of the four vertices in minV/maxV
        const size_t cgindices[4] = { 0, 1, mesh.nx + 1, mesh.nx + 2 };

        assert(alphax.rows() == mesh.nelements);
#pragma omp parallel for schedule(static)
        for (size_t c = 0; c < mesh.nelements; ++c) {
            size_t cx = c % mesh.nx;
            size_t cy = c / mesh.nx;
            size_t cgi = cy * (mesh.nx + 1) + cx; // index of lower-left vertex

            // values of (d/dx phi) in the 4 nodes: lower-left, lower-right, upper-left, upper-right
            const Eigen::Vector<FloatType, 4> vertexvalues = {
                phi(c, 1) - phi(c, 3) - 0.5 * phi(c, 5), phi(c, 1) + phi(c, 3) - 0.5 * phi(c, 5),
                phi(c, 1) - phi(c, 3) + 0.5 * phi(c, 5), phi(c, 1) + phi(c, 3) + 0.5 * phi(c, 5)
            };
            // value of d/dx phi in the midpoint
            const FloatType midvalue = phi(c, 1);

            FloatType al = 1.0; // the limiter
            for (size_t i = 0; i < 4; ++i) {
                FloatType dv = vertexvalues[i] - midvalue; // distance to midpoint
                if (dv > 1.e-8) {
                    assert(_max(cgi + cgindices[i]) >= midvalue);
                    al = std::min(al, std::min(1.0, (_max(cgi + cgindices[i]) - midvalue) / dv));
                }
                if (dv < -1.e-8) {
                    assert(_min(cgi + cgindices[i]) <= midvalue);
                    al = std::min(al, std::min(1.0, (_min(cgi + cgindices[i]) - midvalue) / dv));
                }
                assert(al >= 0);
            }
            alphax(c) = al;
        }
    }
    void computeAlphasY(DGVector<1>& alphay, const DGVector<DG>& phi, const CGVector<1>& _min,
        const CGVector<1>& _max)
    {
        // relative indices of the four vertices in minV/maxV
        const size_t cgindices[4] = { 0, 1, mesh.nx + 1, mesh.nx + 2 };

        assert(alphay.rows() == mesh.nelements);
#pragma omp parallel for schedule(static)
        for (size_t c = 0; c < mesh.nelements; ++c) {
            size_t cx = c % mesh.nx;
            size_t cy = c / mesh.nx;
            size_t cgi = cy * (mesh.nx + 1) + cx; // index of lower-left vertex

            // values of (d/dx phi) in the 4 nodes: lower-left, lower-right, upper-left, upper-right
            const Eigen::Vector<FloatType, 4> vertexvalues = {
                phi(c, 2) - phi(c, 4) - 0.5 * phi(c, 5), phi(c, 2) - phi(c, 4) + 0.5 * phi(c, 5),
                phi(c, 2) + phi(c, 4) - 0.5 * phi(c, 5), phi(c, 2) + phi(c, 4) + 0.5 * phi(c, 5)
            };
            // value of d/dx phi in the midpoint
            const FloatType midvalue = phi(c, 2);

            FloatType al = 1.0; // the limiter
            for (size_t i = 0; i < 4; ++i) {
                FloatType dv = vertexvalues[i] - midvalue; // distance to midpoint
                if (dv > 1.e-8) {
                    assert(_max(cgi + cgindices[i]) >= midvalue);
                    al = std::min(al, std::min(1.0, (_max(cgi + cgindices[i]) - midvalue) / dv));
                }
                if (dv < -1.e-8) {
                    assert(_min(cgi + cgindices[i]) <= midvalue);
                    al = std::min(al, std::min(1.0, (_min(cgi + cgindices[i]) - midvalue) / dv));
                }
                assert(al >= 0);
            }
            alphay(c) = al;
        }
    }

    // performs the vertex-based limiting
    void limit(DGVector<DG>& phi)
    {
        if (DG == 1) // no limiting for dG0
            return;

        // zero order terms & first derivative
        initMinMax(minV, maxV, phi); // get max/min values in vertices
        computeAlphas(alpha, phi, minV, maxV);

        // derivative & second
        if (DG == 6) {
            initMinMax(dxminV, dxmaxV, phi, 1); // get max/min values in vertices
            initMinMax(dyminV, dymaxV, phi, 2); // get max/min values in vertices
            computeAlphasX(alphaX, phi, dxminV, dxmaxV);
            computeAlphasY(alphaY, phi, dyminV, dymaxV);
            for (int i = 0; i < alphaX.rows(); ++i) {
                alphaX(i) = std::min(alphaX(i), alphaY(i));
                alpha(i) = std::max(alpha(i), alphaX(i));
            }
        }

        // limit 2nd order terms
#pragma omp parallel for
        for (int c = 0; c < mesh.nelements; ++c) {
            for (int d = 1; d < 3; ++d)
                phi(c, d) *= alpha(c, 0);
            for (int d = 3; d < DG; ++d)
                phi(c, d) *= alphaX(c, 0);
        }
    }
};

}

#endif
