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
    void InitMinMax(
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
                    const double meanvalue = phi(c, comp);
                    for (size_t j = 0; j < 4; ++j) {
                        _minv(cgi + cgindices[j]) = std::min(_minv(cgi + cgindices[j]), meanvalue);
                        _maxv(cgi + cgindices[j]) = std::max(_maxv(cgi + cgindices[j]), meanvalue);
                    }
                }
            }
    }

    // truncates the averages by min or max value
    void LimitMax(DGVector<DG>& phi, double max) const
    {
#pragma omp parallel for
        for (size_t c = 0; c < mesh.nelements; ++c)
            phi(c, 0) = std::min(max, phi(c, 0));
    }
    // truncates the averages by min or max value
    void LimitMin(DGVector<DG>& phi, double min) const
    {
#pragma omp parallel for
        for (size_t c = 0; c < mesh.nelements; ++c)
            phi(c, 0) = std::max(min, phi(c, 0));
    }

    void ComputeAlphas(DGVector<1>& alpha, const DGVector<DG>& phi, const CGVector<1>& _min,
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
            const Eigen::Vector<double, 4> vertexvalues
                = { phi(c, 0) - 0.5 * phi(c, 1) - 0.5 * phi(c, 2),
                      phi(c, 0) + 0.5 * phi(c, 1) - 0.5 * phi(c, 2),
                      phi(c, 0) - 0.5 * phi(c, 1) + 0.5 * phi(c, 2),
                      phi(c, 0) + 0.5 * phi(c, 1) + 0.5 * phi(c, 2) };

            //	    phi.row(c) * PSILagrange<DG,2>;
            // value of phi in the midpoint
            const double midvalue = phi(c, 0);

            double al = 1.0; // the limiter
            for (size_t i = 0; i < 4; ++i) {
                double dv = vertexvalues[i] - midvalue; // distance to midpoint
                if (dv > 1.e-8) {
                    assert(_max(cgi + cgindices[i]) >= midvalue);
                    al = std::min(al, std::min(1.0, (_max(cgi + cgindices[i]) - midvalue) / dv));
                }
                if (dv < -1.e-8) {
                    // std::cout << std::setprecision(12);
                    // std::cout << minV(cgi + cgindices[i]) << "\t" << midvalue << "\t"
                    //           << vertexvalues[i] << std::endl;
                    assert(_min(cgi + cgindices[i]) <= midvalue);
                    al = std::min(al, std::min(1.0, (_min(cgi + cgindices[i]) - midvalue) / dv));
                }
                assert(al >= 0);
            }
            alpha(c) = al;
        }
    }

    void ComputeAlphasX(DGVector<1>& alphax, const DGVector<DG>& phi, const CGVector<1>& _min,
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
            const Eigen::Vector<double, 4> vertexvalues = { phi(c, 1) - phi(c, 3) - 0.5 * phi(c, 5),
                phi(c, 1) + phi(c, 3) - 0.5 * phi(c, 5), phi(c, 1) - phi(c, 3) + 0.5 * phi(c, 5),
                phi(c, 1) + phi(c, 3) + 0.5 * phi(c, 5) };
            // value of d/dx phi in the midpoint
            const double midvalue = phi(c, 1);

            double al = 1.0; // the limiter
            for (size_t i = 0; i < 4; ++i) {
                double dv = vertexvalues[i] - midvalue; // distance to midpoint
                if (dv > 1.e-8) {
                    assert(_max(cgi + cgindices[i]) >= midvalue);
                    al = std::min(al, std::min(1.0, (_max(cgi + cgindices[i]) - midvalue) / dv));
                }
                if (dv < -1.e-8) {
                    // std::cout << std::setprecision(12);
                    // std::cout << minV(cgi + cgindices[i]) << "\t" << midvalue << "\t"
                    //           << vertexvalues[i] << std::endl;
                    assert(_min(cgi + cgindices[i]) <= midvalue);
                    al = std::min(al, std::min(1.0, (_min(cgi + cgindices[i]) - midvalue) / dv));
                }
                assert(al >= 0);
            }
            alphax(c) = al;
        }
    }
    void ComputeAlphasY(DGVector<1>& alphay, const DGVector<DG>& phi, const CGVector<1>& _min,
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
            const Eigen::Vector<double, 4> vertexvalues = { phi(c, 2) - phi(c, 4) - 0.5 * phi(c, 5),
                phi(c, 2) - phi(c, 4) + 0.5 * phi(c, 5), phi(c, 2) + phi(c, 4) - 0.5 * phi(c, 5),
                phi(c, 2) + phi(c, 4) + 0.5 * phi(c, 5) };
            // value of d/dx phi in the midpoint
            const double midvalue = phi(c, 2);

            double al = 1.0; // the limiter
            for (size_t i = 0; i < 4; ++i) {
                double dv = vertexvalues[i] - midvalue; // distance to midpoint
                if (dv > 1.e-8) {
                    assert(_max(cgi + cgindices[i]) >= midvalue);
                    al = std::min(al, std::min(1.0, (_max(cgi + cgindices[i]) - midvalue) / dv));
                }
                if (dv < -1.e-8) {
                    // std::cout << std::setprecision(12);
                    // std::cout << minV(cgi + cgindices[i]) << "\t" << midvalue << "\t"
                    //           << vertexvalues[i] << std::endl;
                    assert(_min(cgi + cgindices[i]) <= midvalue);
                    al = std::min(al, std::min(1.0, (_min(cgi + cgindices[i]) - midvalue) / dv));
                }
                assert(al >= 0);
            }
            alphay(c) = al;
        }
    }

    // performs the vertex-based limiting
    void Limit(DGVector<DG>& phi)
    {
        if (DG == 1) // no limiting for dG0
            return;

        // zero order terms & first derivative
        InitMinMax(minV, maxV, phi); // get max/min values in vertices
        ComputeAlphas(alpha, phi, minV, maxV);

        // derivative & second
        if (DG == 6) {
            InitMinMax(dxminV, dxmaxV, phi, 1); // get max/min values in vertices
            InitMinMax(dyminV, dymaxV, phi, 2); // get max/min values in vertices
            ComputeAlphasX(alphaX, phi, dxminV, dxmaxV);
            ComputeAlphasY(alphaY, phi, dyminV, dymaxV);
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

    /*
    // performs the vertex-based limiting
    void Limit(DGVector<DG>& phi)
    {
        if (DG == 1) // no limiting for dG0
            return;

        DGVector<3> Dphi(mesh);

        if (DG == 6) // first, higher order limiting
        {
            // limit X-derivative
            Dphi.col(0) = phi.col(1);
            Dphi.col(1) = 2.0 * phi.col(3);
            Dphi.col(2) = phi.col(5);
            InitMinMax(Dphi); // get max/min values in vertices
            ComputeAlphas(alpha, Dphi);

            // limit Y-derivative
            Dphi.col(0) = phi.col(2);
            Dphi.col(1) = phi.col(5);
            Dphi.col(2) = 2.0 * phi.col(4);
            InitMinMax(Dphi); // get max/min values in vertices
            ComputeAlphas(alphaX, Dphi);

            // take minimum alpha for derivative
#pragma omp parallel for
            for (size_t c = 0; c < mesh.nelements; ++c)
                alphaX(c) = std::min(alphaX(c), alpha(c));

            // limit 2nd order terms
#pragma omp parallel for
            for (int c = 0; c < mesh.nelements; ++c) {
                assert(alphaX(c, 0) >= 0.0);
                assert(alphaX(c, 0) <= 1.0);
                if (alphaX(c, 0) < 1) {
                    phi(c, 3) *= alphaX(c, 0) * 0;
                    phi(c, 4) *= alphaX(c, 0) * 0;
                    phi(c, 5) *= alphaX(c, 0) * 0;
                }
            }
        }

        if (DG == 3) // limit first order limiting
        {
            Dphi.col(0) = phi.col(0);
            Dphi.col(1) = phi.col(1);
            Dphi.col(2) = phi.col(2);
            InitMinMax(Dphi); // get max/min values in vertices
            ComputeAlphas(alpha, Dphi);

//             if (DG == 6) // relax?
// #pragma omp parallel for
//                 for (size_t c = 0; c < mesh.nelements; ++c) {
//                     alpha(c) = std::max(alphaX(c), alpha(c));
//                 }
//
//             for (size_t c= 0; c < mesh.nelements; ++c)
//                 if (alpha(c) < 0.99)
//                     std::cout << alpha(c) << std::endl;

            // limit 1st order terms
#pragma omp parallel for
            for (int c = 0; c < mesh.nelements; ++c)
                if (alpha(c, 0) < 1) {
                    phi(c, 1) *= alpha(c, 0);
                    phi(c, 2) *= alpha(c, 0);
                }
        }
    }
    */
};

}

#endif
