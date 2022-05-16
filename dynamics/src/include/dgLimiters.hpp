/*!
 * @file dgLimiters.hpp
 * @date 1 Mar 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __LIMITERS_HPP
#define __LIMITERS_HPP

#include "dgVector.hpp"

namespace Nextsim {

template <int DGdegree>
class Limiter {

    const Mesh& mesh;
    Nextsim::CellVector<DGdegree - 1> alpha;

public:
    Limiter(const Mesh& mesh)
        : mesh(mesh)
        , alpha(mesh)
    {
    }

    void apply_vertex_based_limiter(CellVector<DGdegree>& v)
    {
        // apply limiter
        //#pragma omp parallel for
        for (int ii = 0; ii < v.rows(); ++ii) {
            v(ii, 1) *= alpha(ii, 0);
            v(ii, 2) *= alpha(ii, 0);

            if (DGdegree >= 2) {
                v(ii, 3) *= alpha(ii, 1);
                v(ii, 4) *= alpha(ii, 1);
                v(ii, 5) *= alpha(ii, 1);
            }
        }
    }

    void compute_vertex_based_limiter(const CellVector<1>& v)
    {

        //#pragma omp parallel for
        size_t ii = 0;
        for (size_t iy = 0; iy < mesh.ny; ++iy)
            for (size_t ix = 0; ix < mesh.nx; ++ix, ++ii) {

                // exlude boundaries TODO: implement boundaries
                if ((iy == 0) || (iy == mesh.ny - 1) || (ix == 0) || (ix == mesh.nx - 1))
                    continue;

                // if ((ii % mesh.nx != 0) && ((ii + 1) % mesh.nx != 0) && !(ii <
                // mesh.ny) &&
                //     !(ii > mesh.nx * (mesh.ny - 1))) {

                auto Uc = v(ii, 0);
                // solution values in verticies
                std::array<double, 4> U = { v(ii, 0), v(ii, 0), v(ii, 0), v(ii, 0) };
                U[0] += -0.5 * v(ii, 1);
                U[1] += 0.5 * v(ii, 1);
                U[2] += 0.5 * v(ii, 1);
                U[3] += -0.5 * v(ii, 1);
                U[0] += -0.5 * v(ii, 2);
                U[1] += -0.5 * v(ii, 2);
                U[2] += 0.5 * v(ii, 2);
                U[3] += 0.5 * v(ii, 2);

                // max and min bounds for verticies
                std::array<double, 4> Umin, Umax;

                // loop over verticies taking all elements that contain vertex
                // bottom left, bottom right, top right, top left
                std::array<size_t, 4> shift = { 0, 1, mesh.nx + 1, mesh.nx };
                int j = 0;

                for (auto s : shift) {
                    auto A = { v(ii + s, 0), v(ii - 1 + s, 0), v(ii - mesh.nx - 1 + s, 0),
                        v(ii - mesh.nx + s, 0) };

                    const auto minmax = std::minmax_element(A.begin(), A.end());
                    Umin[j] = *minmax.first;
                    Umax[j] = *minmax.second;
                    j++;
                }

                // compute correction factor alpha
                double alpha_e(1.0);
                for (j = 0; j < 4; j++) {
                    if (std::abs(U[j] - Uc) < 1e-10)
                        alpha_e = std::min(alpha_e, 1.0);
                    else if (U[j] > Uc)
                        alpha_e = std::min(alpha_e, std::min(1.0, (Umax[j] - Uc) / (U[j] - Uc)));
                    else
                        alpha_e = std::min(alpha_e, std::min(1.0, (Umin[j] - Uc) / (U[j] - Uc)));
                }

                alpha(ii, 0) = alpha_e;
            }
    }

    void compute_vertex_based_limiter(const CellVector<2>& v)
    {

        size_t ii = 0;
        for (size_t iy = 0; iy < mesh.ny; ++iy)
            for (size_t ix = 0; ix < mesh.nx; ++ix, ++ii) {
                //    for (size_t ii = 0; ii < v.rows(); ++ii) {

                // exlude boundaries TODO: implement boundaries
                if ((iy == 0) || (iy == mesh.ny - 1) || (ix == 0) || (ix == mesh.nx - 1))
                    continue;

                auto Uc = v(ii, 0);
                auto Ucx = v(ii, 1);
                auto Ucy = v(ii, 2);

                // solution values in verticies
                // Kuzmin2011 eqn 23
                std::array<double, 4> U = { v(ii, 0), v(ii, 0), v(ii, 0), v(ii, 0) };
                U[0] += -0.5 * v(ii, 1);
                U[1] += 0.5 * v(ii, 1);
                U[2] += 0.5 * v(ii, 1);
                U[3] += -0.5 * v(ii, 1);
                U[0] += -0.5 * v(ii, 2);
                U[1] += -0.5 * v(ii, 2);
                U[2] += 0.5 * v(ii, 2);
                U[3] += 0.5 * v(ii, 2);

                // Kuzmin2011 eqns 21 and 22
                std::array<double, 4> Ux = { v(ii, 1), v(ii, 1), v(ii, 1), v(ii, 1) };
                Ux[0] += -0.5 * v(ii, 3);
                Ux[1] += 0.5 * v(ii, 3);
                Ux[2] += 0.5 * v(ii, 3);
                Ux[3] += -0.5 * v(ii, 3);
                Ux[0] += -0.5 * v(ii, 5);
                Ux[1] += -0.5 * v(ii, 5);
                Ux[2] += 0.5 * v(ii, 5);
                Ux[3] += 0.5 * v(ii, 5);
                std::array<double, 4> Uy = { v(ii, 2), v(ii, 2), v(ii, 2), v(ii, 2) };
                Uy[0] += -0.5 * v(ii, 5);
                Uy[1] += 0.5 * v(ii, 5);
                Uy[2] += 0.5 * v(ii, 5);
                Uy[3] += -0.5 * v(ii, 5);
                Uy[0] += -0.5 * v(ii, 4);
                Uy[1] += -0.5 * v(ii, 4);
                Uy[2] += 0.5 * v(ii, 4);
                Uy[3] += 0.5 * v(ii, 4);

                // max and min bounds for verticies
                std::array<double, 4> Umin, Umax, Uminx, Umaxx, Uminy, Umaxy;

                // loop over verticies taking all elements that contain vertex
                // bottom left, bottom right, top right, top left
                std::array<size_t, 4> shift = { 0, 1, mesh.nx + 1, mesh.nx };
                int j = 0;

                for (auto s : shift) {
                    auto A = { v(ii + s, 0), v(ii - 1 + s, 0), v(ii - mesh.nx - 1 + s, 0),
                        v(ii - mesh.nx + s, 0) };
                    auto Ax = { v(ii + s, 1), v(ii - 1 + s, 1), v(ii - mesh.nx - 1 + s, 1),
                        v(ii - mesh.nx + s, 1) };
                    auto Ay = { v(ii + s, 2), v(ii - 1 + s, 2), v(ii - mesh.nx - 1 + s, 2),
                        v(ii - mesh.nx + s, 2) };

                    const auto minnax = std::minmax_element(A.begin(), A.end());
                    const auto minxmaxx = std::minmax_element(Ax.begin(), Ax.end());
                    const auto minymaxy = std::minmax_element(Ay.begin(), Ay.end());
                    Umin[j] = *minnax.first;
                    Umax[j] = *minnax.second;
                    Uminx[j] = *minxmaxx.first;
                    Umaxx[j] = *minxmaxx.second;
                    Uminy[j] = *minymaxy.first;
                    Umaxy[j] = *minymaxy.second;

                    j++;
                }

                // compute correction factor alpha
                double alpha2x(1.0);
                for (j = 0; j < 4; j++) {
                    if (std::abs(Ux[j] - Ucx) < 1e-10)
                        alpha2x = std::min(alpha2x, 1.0);
                    else if (Ux[j] > Ucx)
                        alpha2x
                            = std::min(alpha2x, std::min(1.0, (Umaxx[j] - Ucx) / (Ux[j] - Ucx)));
                    else
                        alpha2x
                            = std::min(alpha2x, std::min(1.0, (Uminx[j] - Ucx) / (Ux[j] - Ucx)));
                }
                double alpha2y(1.0);
                for (j = 0; j < 4; j++) {
                    if (std::abs(Uy[j] - Ucy) < 1e-10)
                        alpha2y = std::min(alpha2y, 1.0);
                    else if (Uy[j] > Ucy)
                        alpha2y
                            = std::min(alpha2y, std::min(1.0, (Umaxy[j] - Ucy) / (Uy[j] - Ucy)));
                    else
                        alpha2y
                            = std::min(alpha2y, std::min(1.0, (Uminy[j] - Ucy) / (Uy[j] - Ucy)));
                }

                alpha(ii, 1) = std::min(alpha2x, alpha2y);

                double alpha_e(1.0);
                for (j = 0; j < 4; j++) {
                    if (std::abs(U[j] - Uc) < 1e-10)
                        alpha_e = std::min(alpha_e, 1.0);
                    else if (U[j] > Uc)
                        alpha_e = std::min(alpha_e, std::min(1.0, (Umax[j] - Uc) / (U[j] - Uc)));
                    else
                        alpha_e = std::min(alpha_e, std::min(1.0, (Umin[j] - Uc) / (U[j] - Uc)));
                }

                alpha(ii, 0) = std::max(alpha_e, alpha(ii, 1));
            }
    }
};
}

#endif /* __LIMITERS_HPP */
