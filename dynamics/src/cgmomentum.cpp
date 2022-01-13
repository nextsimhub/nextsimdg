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
            Eigen::Vector<double, 4> cg_local = {
                cg(cgi), cg(cgi + 1),
                cg(cgi + cgshift), cg(cgi + 1 + cgshift)
            };

            dg.row(dgi) = CG1_to_DG2 * cg_local;
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
            Eigen::Vector<double, 4> vx_local = {
                vx(cgi), vx(cgi + 1),
                vx(cgi + cgshift), vx(cgi + 1 + cgshift)
            };
            Eigen::Vector<double, 4> vy_local = {
                vy(cgi), vy(cgi + 1),
                vy(cgi + cgshift), vy(cgi + 1 + cgshift)
            };

            E11.row(dgi) = CG1_to_DG1_dX * vx_local / mesh.hx;
            E22.row(dgi) = CG1_to_DG1_dY * vy_local / mesh.hy;
            E12.row(dgi) = 0.5 * CG1_to_DG1_dX * vy_local / mesh.hx + 0.5 * CG1_to_DG1_dY * vx_local / mesh.hy;
        }
    }
}

//! += scale * (S, nabla Phi)  with Phi=(PhiX, PhiY) CG2 - test functions
template <>
void CGMomentum::AddStressTensor(const Mesh& mesh, double scale,
    CGVector<2>& tx, CGVector<2>& ty,
    const CellVector<1>& S11, const CellVector<1>& S12, const CellVector<1>& S22) const
{
    const size_t CGROW = 2 * mesh.nx + 1;

    //    //! For parallelization and memory access control use simple checkerboard pattern
    size_t i = 0;
    for (size_t iy = 0; iy < mesh.ny; ++iy)
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++i) {

            //     for (size_t px = 0; px < 2; ++px)
            //         for (size_t py = 0; py < 2; ++py) {

            // #pragma omp parallel for
            //            for (size_t i = 0; i < mesh.n; ++i) {
            // if (ix % 2 == px) //!< Only proceed if pattern matches.
            //     continue;
            // if (iy % 2 == py)
            //     continue;

            size_t cg_i = 2 * CGROW * iy + 2 * ix; //!< Lower left index of element in CG2-Vector

            Eigen::Matrix<double, 9, 1> lup1 = scale * (DG1_CG2_dX * S11.row(i).transpose() / mesh.hx + DG1_CG2_dY * S12.row(i).transpose() / mesh.hy);
            tx(cg_i + 0) += lup1(0);
            tx(cg_i + 1) += lup1(1);
            tx(cg_i + 2) += lup1(2);
            tx(cg_i + 0 + CGROW) += lup1(3);
            tx(cg_i + 1 + CGROW) += lup1(4);
            tx(cg_i + 2 + CGROW) += lup1(5);
            tx(cg_i + 0 + CGROW * 2) += lup1(6);
            tx(cg_i + 1 + CGROW * 2) += lup1(7);
            tx(cg_i + 2 + CGROW * 2) += lup1(8);

            Eigen::Matrix<double, 9, 1> lup2 = scale * (DG1_CG2_dX * S12.row(i).transpose() / mesh.hx + DG1_CG2_dY * S22.row(i).transpose() / mesh.hy);
            ty(cg_i + 0) += lup2(0);
            ty(cg_i + 1) += lup2(1);
            ty(cg_i + 2) += lup2(2);
            ty(cg_i + 0 + CGROW) += lup2(3);
            ty(cg_i + 1 + CGROW) += lup2(4);
            ty(cg_i + 2 + CGROW) += lup2(5);
            ty(cg_i + 0 + CGROW * 2) += lup2(6);
            ty(cg_i + 1 + CGROW * 2) += lup2(7);
            ty(cg_i + 2 + CGROW * 2) += lup2(8);
        }
}
template <>
void CGMomentum::AddStressTensor(const Mesh& mesh, double scale,
    CGVector<1>& tx, CGVector<1>& ty,
    const CellVector<1>& S11, const CellVector<1>& S12, const CellVector<1>& S22) const
{
    const size_t CGROW = mesh.nx + 1;

    //    //! For parallelization and memory access control use simple checkerboard pattern
    size_t i = 0;
    for (size_t iy = 0; iy < mesh.ny; ++iy)
        for (size_t ix = 0; ix < mesh.nx; ++ix, ++i) {

            //     for (size_t px = 0; px < 2; ++px)
            //         for (size_t py = 0; py < 2; ++py) {

            // #pragma omp parallel for
            //            for (size_t i = 0; i < mesh.n; ++i) {
            // if (ix % 2 == px) //!< Only proceed if pattern matches.
            //     continue;
            // if (iy % 2 == py)
            //     continue;

            size_t cg_i = CGROW * iy + ix; //!< Lower left index of element in CG2-Vector

            Eigen::Matrix<double, 4, 1> lup1 = scale * (DG1_CG1_dX * S11.row(i).transpose() / mesh.hx + DG1_CG1_dY * S12.row(i).transpose() / mesh.hy);
            tx(cg_i + 0) += lup1(0);
            tx(cg_i + 1) += lup1(1);
            tx(cg_i + 0 + CGROW) += lup1(2);
            tx(cg_i + 1 + CGROW) += lup1(3);

            Eigen::Matrix<double, 4, 1> lup2 = scale * (DG1_CG1_dX * S12.row(i).transpose() / mesh.hx + DG1_CG1_dY * S22.row(i).transpose() / mesh.hy);
            ty(cg_i + 0) += lup2(0);
            ty(cg_i + 1) += lup2(1);
            ty(cg_i + 0 + CGROW) += lup2(2);
            ty(cg_i + 1 + CGROW) += lup2(3);
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

}
