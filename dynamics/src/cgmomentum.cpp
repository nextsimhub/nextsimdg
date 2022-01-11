#include "cgmomentum.hpp"
#include "codegeneration_cg_to_dg.hpp"

namespace Nextsim {

template <>
void CGMomentum::ProjectCGToDG(const Mesh& mesh, CellVector<1>& dg, const CGVector<2>& cg)
{
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = 2 * mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 2) {
            Eigen::Vector<double, 9> cg_local = {
                cg(cgi), cg(cgi + 1), cg(cgi + 2),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift),
                cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift)
            };
            dg.row(dgi) = CG2_to_DG1 * cg_local;
        }
    }
}
template <>
void CGMomentum::ProjectCGToDG(const Mesh& mesh, CellVector<2>& dg, const CGVector<2>& cg)
{
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = 2 * mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 2) {
            Eigen::Vector<double, 9> cg_local = {
                cg(cgi), cg(cgi + 1), cg(cgi + 2),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift),
                cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift)
            };

            dg.row(dgi) = CG2_to_DG2 * cg_local;
        }
    }
}

/*!
   * projects the symmatric gradient of the continuous CG2 velocity into a dg vector
   */
void CGMomentum::ProjectCG2VelocityToDG1Strain(const Mesh& mesh,
    CellVector<1>& E11, CellVector<1>& E12, CellVector<1>& E22,
    const CGVector<2>& vx, const CGVector<2>& vy)
{
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == vx.rows());
    assert(static_cast<long int>((2 * mesh.nx + 1) * (2 * mesh.ny + 1)) == vy.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E11.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E12.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E22.rows());

    const int cgshift = 2 * mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 2) {
            Eigen::Vector<double, 9> vx_local = {
                vx(cgi), vx(cgi + 1), vx(cgi + 2),
                vx(cgi + cgshift), vx(cgi + 1 + cgshift), vx(cgi + 2 + cgshift),
                vx(cgi + 2 * cgshift), vx(cgi + 1 + 2 * cgshift), vx(cgi + 2 + 2 * cgshift)
            };
            Eigen::Vector<double, 9> vy_local = {
                vy(cgi), vy(cgi + 1), vy(cgi + 2),
                vy(cgi + cgshift), vy(cgi + 1 + cgshift), vy(cgi + 2 + cgshift),
                vy(cgi + 2 * cgshift), vy(cgi + 1 + 2 * cgshift), vy(cgi + 2 + 2 * cgshift)
            };

            E11.row(dgi) = CG2_to_DG1_dX * vx_local / mesh.hx;
            E22.row(dgi) = CG2_to_DG1_dY * vy_local / mesh.hy;
            E12.row(dgi) = 0.5 * CG2_to_DG1_dX * vy_local / mesh.hx + 0.5 * CG2_to_DG1_dY * vx_local / mesh.hy;
        }
    }
}

}
