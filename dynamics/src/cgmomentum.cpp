#include "cgmomentum.hpp"

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
            Eigen::Matrix<double, 9, 1> cg_local;
            cg_local << cg(cgi), cg(cgi + 1), cg(cgi + 2),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift),
                cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift);
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
            Eigen::Matrix<double, 9, 1> cg_local;
            cg_local << cg(cgi), cg(cgi + 1), cg(cgi + 2),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift), cg(cgi + 2 + cgshift),
                cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift), cg(cgi + 2 + 2 * cgshift);

            dg.row(dgi) = CG2_to_DG2 * cg_local;
        }
    }
}
template <>
void CGMomentum::ProjectCGToDG(const Mesh& mesh, CellVector<2>& dg, const CGVector<1>& cg)
{
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 1) {
            Eigen::Matrix<double, 4, 1> cg_local;
            cg_local << cg(cgi), cg(cgi + 1),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift);

            dg.row(dgi) = CG1_to_DG2 * cg_local;
        }
    }
}
template <>
void CGMomentum::ProjectCGToDG(const Mesh& mesh, CellVector<1>& dg, const CGVector<1>& cg)
{
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 1) {
            Eigen::Matrix<double, 4, 1> cg_local;
            cg_local << cg(cgi), cg(cgi + 1),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift);

            dg.row(dgi) = CG1_to_DG1 * cg_local;
        }
    }
}
template <>
void CGMomentum::ProjectCGToDG(const Mesh& mesh, CellVector<0>& dg, const CGVector<1>& cg)
{
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == dg.rows());

    const int cgshift = mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 1) {
            Eigen::Matrix<double, 4, 1> cg_local;
            cg_local << cg(cgi), cg(cgi + 1),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift);

            dg.row(dgi) = CG1_to_DG0 * cg_local;
        }
    }
}

/*!
   * projects the symmatric gradient of the continuous CG2 velocity into a dg vector
   */
template <>
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
            Eigen::Matrix<double, 9, 1> vx_local;
            vx_local << vx(cgi), vx(cgi + 1), vx(cgi + 2),
                vx(cgi + cgshift), vx(cgi + 1 + cgshift), vx(cgi + 2 + cgshift),
                vx(cgi + 2 * cgshift), vx(cgi + 1 + 2 * cgshift), vx(cgi + 2 + 2 * cgshift);
            Eigen::Matrix<double, 9, 1> vy_local;
            vy_local << vy(cgi), vy(cgi + 1), vy(cgi + 2),
                vy(cgi + cgshift), vy(cgi + 1 + cgshift), vy(cgi + 2 + cgshift),
                vy(cgi + 2 * cgshift), vy(cgi + 1 + 2 * cgshift), vy(cgi + 2 + 2 * cgshift);

            E11.row(dgi) = CG2_to_DG1_dX * vx_local / mesh.hx;
            E22.row(dgi) = CG2_to_DG1_dY * vy_local / mesh.hy;
            E12.row(dgi) = 0.5 * CG2_to_DG1_dX * vy_local / mesh.hx + 0.5 * CG2_to_DG1_dY * vx_local / mesh.hy;
        }
    }
}
template <>
void CGMomentum::ProjectCG2VelocityToDG1Strain(const Mesh& mesh,
    CellVector<1>& E11, CellVector<1>& E12, CellVector<1>& E22,
    const CGVector<1>& vx, const CGVector<1>& vy)
{
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == vx.rows());
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == vy.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E11.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E12.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E22.rows());

    const int cgshift = mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 1) {
            Eigen::Matrix<double, 4, 1> vx_local = {
                vx(cgi), vx(cgi + 1),
                vx(cgi + cgshift), vx(cgi + 1 + cgshift)
            };
            Eigen::Matrix<double, 4, 1> vy_local = {
                vy(cgi), vy(cgi + 1),
                vy(cgi + cgshift), vy(cgi + 1 + cgshift)
            };

            E11.row(dgi) = CG1_to_DG1_dX * vx_local / mesh.hx;
            E22.row(dgi) = CG1_to_DG1_dY * vy_local / mesh.hy;
            E12.row(dgi) = 0.5 * CG1_to_DG1_dX * vy_local / mesh.hx + 0.5 * CG1_to_DG1_dY * vx_local / mesh.hy;
        }
    }
}
template <>
void CGMomentum::ProjectCG2VelocityToDG1Strain(const Mesh& mesh,
    CellVector<0>& E11, CellVector<0>& E12, CellVector<0>& E22,
    const CGVector<1>& vx, const CGVector<1>& vy)
{
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == vx.rows());
    assert(static_cast<long int>((mesh.nx + 1) * (mesh.ny + 1)) == vy.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E11.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E12.rows());
    assert(static_cast<long int>(mesh.nx * mesh.ny) == E22.rows());

    const int cgshift = mesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < mesh.ny; ++row) {
        int dgi = mesh.nx * row; //!< Index of dg vector
        int cgi = cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < mesh.nx; ++col, ++dgi, cgi += 1) {
            Eigen::Matrix<double, 4, 1> vx_local = {
                vx(cgi), vx(cgi + 1),
                vx(cgi + cgshift), vx(cgi + 1 + cgshift)
            };
            Eigen::Matrix<double, 4, 1> vy_local = {
                vy(cgi), vy(cgi + 1),
                vy(cgi + cgshift), vy(cgi + 1 + cgshift)
            };

            E11.row(dgi) = CG1_to_DG0_dX * vx_local / mesh.hx;
            E22.row(dgi) = CG1_to_DG0_dY * vy_local / mesh.hy;
            E12.row(dgi) = 0.5 * CG1_to_DG0_dX * vy_local / mesh.hx + 0.5 * CG1_to_DG0_dY * vx_local / mesh.hy;
        }
    }
}

template <int CG, int DG>
void CGMomentum::AddStressTensor(const Mesh& mesh, const double scale,
    CGVector<CG>& tx, CGVector<CG>& ty,
    const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const
{
    // parallelization in tripes
    for (size_t p = 0; p < 2; ++p)
#pragma omp parallel for schedule(static)
        for (size_t cy = 0; cy < mesh.ny; ++cy) //!< loop over all cells of the mesh
        {
            if (cy % 2 == p) {
                size_t c = mesh.nx * cy;
                for (size_t cx = 0; cx < mesh.nx; ++cx, ++c) //!< loop over all cells of the mesh
                    AddStressTensorCell(mesh, scale, c, cx, cy, tx, ty, S11, S12, S22);
            }
        }
}

//! Sets the vector to zero along the boundary
template <>
void CGMomentum::DirichletZero(const Mesh& mesh, CGVector<1>& v) const
{

    size_t upperleftindex = (mesh.nx + 1) * mesh.ny;
    for (size_t i = 0; i < mesh.nx + 1; ++i) {
        v(i, 0) = 0.0;
        v(upperleftindex + i, 0) = 0.0;
    }
    size_t indecesperrow = mesh.nx + 1;
    size_t lowerrightindex = mesh.nx;
    for (size_t i = 0; i < mesh.ny + 1; ++i) {
        v(indecesperrow * i, 0) = 0.0;
        v(lowerrightindex + indecesperrow * i, 0) = 0.0;
    }
}
template <>
void CGMomentum::DirichletZero(const Mesh& mesh, CGVector<2>& v) const
{

    size_t upperleftindex = (2 * mesh.nx + 1) * 2 * mesh.ny;
    for (size_t i = 0; i < 2 * mesh.nx + 1; ++i) {
        v(i, 0) = 0.0;
        v(upperleftindex + i, 0) = 0.0;
    }
    size_t indecesperrow = 2 * mesh.nx + 1;
    size_t lowerleftindex = 2 * mesh.nx;
    for (size_t i = 0; i < 2 * mesh.ny + 1; ++i) {
        v(indecesperrow * i, 0) = 0.0;
        v(lowerleftindex + indecesperrow * i, 0) = 0.0;
    }
}

template <int CG, int DG>
void CGMomentum::InterpolateDGToCG(const Mesh& mesh, CGVector<CG>& cg_A, const CellVector<DG>& A) const
{
    cg_A.zero();

    // parallelization by running over stripes
    for (size_t p = 0; p < 2; ++p) {
#pragma omp parallel for
        for (size_t cy = 0; cy < mesh.ny; ++cy) {
            if (cy % 2 == p)
                continue;

            size_t c = cy * mesh.nx;

            for (size_t cx = 0; cx < mesh.nx; ++cx, ++c)
                InterpolateDGToCGCell(mesh, c, cx, cy, cg_A, A);
        }
    }

    // InterpolateDGToCGCell adds to the CG-Dofs with correct weighting. Then we adjust the boundary
    InterpolateDGToCGBoundary(mesh, cg_A);
}

// --------------------------------------------------

template void CGMomentum::AddStressTensor(const Mesh& mesh, const double scale,
    CGVector<1>& tx, CGVector<1>& ty,
    const CellVector<0>& S11, const CellVector<0>& S12, const CellVector<0>& S22) const;

template void CGMomentum::AddStressTensor(const Mesh& mesh, const double scale,
    CGVector<1>& tx, CGVector<1>& ty,
    const CellVector<1>& S11, const CellVector<1>& S12, const CellVector<1>& S22) const;

template void CGMomentum::AddStressTensor(const Mesh& mesh, const double scale,
    CGVector<2>& tx, CGVector<2>& ty,
    const CellVector<1>& S11, const CellVector<1>& S12, const CellVector<1>& S22) const;

template void CGMomentum::InterpolateDGToCG(const Mesh& mesh, CGVector<1>& cg_A, const CellVector<0>& A) const;
template void CGMomentum::InterpolateDGToCG(const Mesh& mesh, CGVector<2>& cg_A, const CellVector<0>& A) const;

}
