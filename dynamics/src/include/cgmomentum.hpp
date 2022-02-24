/*----------------------------   cgmomentum.hpp     ---------------------------*/
#ifndef __cgmomentum_HPP
#define __cgmomentum_HPP
/*----------------------------   cgmomentum.hpp     ---------------------------*/

#include "cgvector.hpp"
#include "codegeneration_cg_to_dg.hpp"
#include "dgvector.hpp"

namespace Nextsim {

class CGMomentum {
private:
public:
    /*!
     * The following functions take care of the interpolation and projection
     * between CG and DG functions
     */
    //! Projects a CG function to a DG function
    template <int CG, int DG>
    void ProjectCGToDG(const Mesh& mesh, CellVector<DG>& dg, const CGVector<CG>& cg);

    //! Projects the symmetric gradient of the CG2 velocity into the DG1 space
    template <int CG, int DG>
    void ProjectCG2VelocityToDG1Strain(const Mesh& mesh,
        CellVector<DG>& E11, CellVector<DG>& E12, CellVector<DG>& E22,
        const CGVector<CG>& vx, const CGVector<CG>& vy);

    /*!
     * Evaluates (S, nabla phi) and adds it to tx/ty - Vector
     */
    template <int CG, int DG>
    void AddStressTensor(const Mesh& mesh, const double scale,
        CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    template <int CG, int DG>
    void AddStressTensorCell(const Mesh& mesh, const double scale,
        const size_t c, const size_t cx, const size_t cy,
        CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    void AddStressTensorCell(const Mesh& mesh, const double scale,
        const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& tx, CGVector<1>& ty,
        const CellVector<0>& S11, const CellVector<0>& S12, const CellVector<0>& S22) const
    {
        const size_t CGROW = mesh.nx + 1;
        const size_t cg_i = CGROW * cy + cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 4, 1> lup1 = scale * (DG0_CG1_dX * S11.row(c).transpose() / mesh.hx + DG0_CG1_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 0 + CGROW) += lup1(2);
        tx(cg_i + 1 + CGROW) += lup1(3);

        Eigen::Matrix<double, 4, 1> lup2 = scale * (DG0_CG1_dX * S12.row(c).transpose() / mesh.hx + DG0_CG1_dY * S22.row(c).transpose() / mesh.hy);
        ty(cg_i + 0) += lup2(0);
        ty(cg_i + 1) += lup2(1);
        ty(cg_i + 0 + CGROW) += lup2(2);
        ty(cg_i + 1 + CGROW) += lup2(3);
    }
    void AddStressTensorCell(const Mesh& mesh, const double scale,
        const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& tx, CGVector<1>& ty,
        const CellVector<1>& S11, const CellVector<1>& S12, const CellVector<1>& S22) const
    {
        const size_t CGROW = mesh.nx + 1;
        const size_t cg_i = CGROW * cy + cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 4, 1> lup1 = scale * (DG1_CG1_dX * S11.row(c).transpose() / mesh.hx + DG1_CG1_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 0 + CGROW) += lup1(2);
        tx(cg_i + 1 + CGROW) += lup1(3);

        Eigen::Matrix<double, 4, 1> lup2 = scale * (DG1_CG1_dX * S12.row(c).transpose() / mesh.hx + DG1_CG1_dY * S22.row(c).transpose() / mesh.hy);
        ty(cg_i + 0) += lup2(0);
        ty(cg_i + 1) += lup2(1);
        ty(cg_i + 0 + CGROW) += lup2(2);
        ty(cg_i + 1 + CGROW) += lup2(3);
    }
    void AddStressTensorCell(const Mesh& mesh, const double scale,
        const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& tx, CGVector<2>& ty,
        const CellVector<1>& S11, const CellVector<1>& S12, const CellVector<1>& S22) const
    {
        const size_t CGROW = 2 * mesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 9, 1> lup1 = scale * (DG1_CG2_dX * S11.row(c).transpose() / mesh.hx + DG1_CG2_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 2) += lup1(2);
        tx(cg_i + 0 + CGROW) += lup1(3);
        tx(cg_i + 1 + CGROW) += lup1(4);
        tx(cg_i + 2 + CGROW) += lup1(5);
        tx(cg_i + 0 + CGROW * 2) += lup1(6);
        tx(cg_i + 1 + CGROW * 2) += lup1(7);
        tx(cg_i + 2 + CGROW * 2) += lup1(8);

        Eigen::Matrix<double, 9, 1> lup2 = scale * (DG1_CG2_dX * S12.row(c).transpose() / mesh.hx + DG1_CG2_dY * S22.row(c).transpose() / mesh.hy);
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

    //! Sets the vector to zero along the boundary
    template <int CG>
    void DirichletZero(const Mesh& mesh, CGVector<CG>& v) const;

    //! Interpolates a DG-Vector to a CG-Vector
    template <int CG, int DG>
    void InterpolateDGToCG(const Mesh& mesh, CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    template <int CG, int DG>
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy, CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy, CGVector<1>& cg_A, const CellVector<0>& A) const
    {
        const size_t CGDofsPerRow = mesh.nx + 1;
        const size_t cgi = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * A(c);
        cg_A(cgi + 1) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow + 1) += 0.25 * A(c);
    }
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy, CGVector<1>& cg_A, const CellVector<1>& A) const
    {
        const size_t CGDofsPerRow = mesh.nx + 1;
        const size_t cgi = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow + 1) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    }
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy, CGVector<2>& cg_A, const CellVector<0>& A) const
    {
        const size_t CGDofsPerRow = 2 * mesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * A(c);
        cg_A(cgi + 1) += 0.5 * A(c);
        cg_A(cgi + 2) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow) += 0.5 * A(c);
        cg_A(cgi + CGDofsPerRow + 1) += A(c);
        cg_A(cgi + CGDofsPerRow + 2) += 0.5 * A(c);
        cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * A(c);
        cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * A(c);
        cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * A(c);
    }

    //! Adjusts the interpolation on the boundary
    template <int CG>
    void InterpolateDGToCGBoundary(const Mesh& mesh, CGVector<CG>& cg_A) const;
    void InterpolateDGToCGBoundary(const Mesh& mesh, CGVector<1>& cg_A) const
    {
        const size_t CGDofsPerRow = mesh.nx + 1;
        const size_t UpperLeftIndex = CGDofsPerRow * mesh.ny;

#pragma omp parallel for
        for (size_t i = 0; i < mesh.nx + 1; ++i) {
            cg_A(i) *= 2.0;
            cg_A(UpperLeftIndex + i) *= 2.0;
        }
#pragma omp parallel for
        for (size_t i = 0; i < mesh.ny + 1; ++i) {
            cg_A(i * CGDofsPerRow) *= 2.0;
            cg_A(i * CGDofsPerRow + mesh.nx) *= 2.0;
        }
    }
    void InterpolateDGToCGBoundary(const Mesh& mesh, CGVector<2>& cg_A) const
    {
        const size_t CGDofsPerRow = 2 * mesh.nx + 1;
        const size_t UpperLeftIndex = 2 * CGDofsPerRow * mesh.ny;

#pragma omp parallel for
        for (size_t i = 0; i < 2 * mesh.nx + 1; ++i) {
            cg_A(i) *= 2.0;
            cg_A(UpperLeftIndex + i) *= 2.0;
        }
#pragma omp parallel for
        for (size_t i = 0; i < 2 * mesh.ny + 1; ++i) {
            cg_A(i * CGDofsPerRow) *= 2.0;
            cg_A(i * CGDofsPerRow + 2 * mesh.nx) *= 2.0;
        }
    }
};

}

/*----------------------------   cgmomentum.hpp     ---------------------------*/
/* end of #ifndef __cgmomentum_HPP */
#endif
/*----------------------------   cgmomentum.hpp     ---------------------------*/
