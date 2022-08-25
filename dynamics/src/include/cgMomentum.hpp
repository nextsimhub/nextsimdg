/*!
 * @file Mesh.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __CGMOMENTUM_HPP
#define __CGMOMENTUM_HPP

#include "ParametricTools.hpp"
#include "cgVector.hpp"
#include "codeGenerationCGToDG.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

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
    void ProjectCG2VelocityToDG1Strain(const Mesh& mesh, CellVector<DG>& E11, CellVector<DG>& E12,
        CellVector<DG>& E22, const CGVector<CG>& vx, const CGVector<CG>& vy);

    //! Projects the symmetric gradient of the CG2 velocity into the DG1 space
    template <int CG, int DG>
    void ProjectCG2VelocityToDG1Strain(const SasipMesh& smesh, CellVector<DG>& E11, CellVector<DG>& E12,
        CellVector<DG>& E22, const CGVector<CG>& vx, const CGVector<CG>& vy);
    template <int CG, int DG>
    void ProjectCG2VelocityToDG1Strain(const ParametricTransformation<CG, DG>& ptrans,
        const SasipMesh& smesh, CellVector<DG>& E11, CellVector<DG>& E12,
        CellVector<DG>& E22, const CGVector<CG>& vx, const CGVector<CG>& vy);

    /*!
     * Evaluates (S, nabla phi) and adds it to tx/ty - Vector
     */
    template <int CG, int DG>
    void AddStressTensor(const Mesh& mesh, const double scale, CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    template <int CG, int DG>
    void AddStressTensor(const SasipMesh& smesh, const double scale, CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    template <int CG, int DG>
    void AddStressTensor(const ParametricTransformation<CG, DG>& ptrans,
        const SasipMesh& smesh, const double scale, CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    template <int CG, int DG>
    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<CG>& tx, CGVector<CG>& ty, const CellVector<DG>& S11,
        const CellVector<DG>& S12, const CellVector<DG>& S22) const;
    template <int CG, int DG>
    void AddStressTensorCell(const SasipMesh& smesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<CG>& tx, CGVector<CG>& ty, const CellVector<DG>& S11,
        const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    template <int CG, int DG>
    void AddStressTensorCell(const ParametricTransformation<CG, DG>& ptrans,
        const SasipMesh& smesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<CG>& tx, CGVector<CG>& ty, const CellVector<DG>& S11,
        const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<1>& tx, CGVector<1>& ty, const CellVector<1>& S11,
        const CellVector<1>& S12, const CellVector<1>& S22) const
    {
        const size_t CGROW = mesh.nx + 1;
        const size_t cg_i = CGROW * cy + cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 4, 1> lup1 = scale
            * (DG1_CG1_dX * S11.row(c).transpose() / mesh.hx
                + DG1_CG1_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 0 + CGROW) += lup1(2);
        tx(cg_i + 1 + CGROW) += lup1(3);

        Eigen::Matrix<double, 4, 1> lup2 = scale
            * (DG1_CG1_dX * S12.row(c).transpose() / mesh.hx
                + DG1_CG1_dY * S22.row(c).transpose() / mesh.hy);
        ty(cg_i + 0) += lup2(0);
        ty(cg_i + 1) += lup2(1);
        ty(cg_i + 0 + CGROW) += lup2(2);
        ty(cg_i + 1 + CGROW) += lup2(3);
    }
    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<1>& tx, CGVector<1>& ty, const CellVector<3>& S11,
        const CellVector<3>& S12, const CellVector<3>& S22) const
    {
        const size_t CGROW = mesh.nx + 1;
        const size_t cg_i = CGROW * cy + cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 4, 1> lup1 = scale
            * (DG3_CG1_dX * S11.row(c).transpose() / mesh.hx
                + DG3_CG1_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 0 + CGROW) += lup1(2);
        tx(cg_i + 1 + CGROW) += lup1(3);

        Eigen::Matrix<double, 4, 1> lup2 = scale
            * (DG3_CG1_dX * S12.row(c).transpose() / mesh.hx
                + DG3_CG1_dY * S22.row(c).transpose() / mesh.hy);
        ty(cg_i + 0) += lup2(0);
        ty(cg_i + 1) += lup2(1);
        ty(cg_i + 0 + CGROW) += lup2(2);
        ty(cg_i + 1 + CGROW) += lup2(3);
    }
    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<1>& tx, CGVector<1>& ty, const CellVector<6>& S11,
        const CellVector<6>& S12, const CellVector<6>& S22) const
    {
        const size_t CGROW = mesh.nx + 1;
        const size_t cg_i = CGROW * cy + cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 4, 1> lup1 = scale
            * (DG6_CG1_dX * S11.row(c).transpose() / mesh.hx
                + DG6_CG1_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 0 + CGROW) += lup1(2);
        tx(cg_i + 1 + CGROW) += lup1(3);

        Eigen::Matrix<double, 4, 1> lup2 = scale
            * (DG6_CG1_dX * S12.row(c).transpose() / mesh.hx
                + DG6_CG1_dY * S22.row(c).transpose() / mesh.hy);
        ty(cg_i + 0) += lup2(0);
        ty(cg_i + 1) += lup2(1);
        ty(cg_i + 0 + CGROW) += lup2(2);
        ty(cg_i + 1 + CGROW) += lup2(3);
    }
    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<2>& tx, CGVector<2>& ty, const CellVector<3>& S11,
        const CellVector<3>& S12, const CellVector<3>& S22) const
    {
        const size_t CGROW = 2 * mesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 9, 1> lup1 = scale
            * (DG3_CG2_dX * S11.row(c).transpose() / mesh.hx
                + DG3_CG2_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 2) += lup1(2);
        tx(cg_i + 0 + CGROW) += lup1(3);
        tx(cg_i + 1 + CGROW) += lup1(4);
        tx(cg_i + 2 + CGROW) += lup1(5);
        tx(cg_i + 0 + CGROW * 2) += lup1(6);
        tx(cg_i + 1 + CGROW * 2) += lup1(7);
        tx(cg_i + 2 + CGROW * 2) += lup1(8);

        Eigen::Matrix<double, 9, 1> lup2 = scale
            * (DG3_CG2_dX * S12.row(c).transpose() / mesh.hx
                + DG3_CG2_dY * S22.row(c).transpose() / mesh.hy);
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

    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<2>& tx, CGVector<2>& ty, const CellVector<6>& S11,
        const CellVector<6>& S12, const CellVector<6>& S22) const
    {
        const size_t CGROW = 2 * mesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 9, 1> lup1 = scale
            * (DG6_CG2_dX * S11.row(c).transpose() / mesh.hx
                + DG6_CG2_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 2) += lup1(2);
        tx(cg_i + 0 + CGROW) += lup1(3);
        tx(cg_i + 1 + CGROW) += lup1(4);
        tx(cg_i + 2 + CGROW) += lup1(5);
        tx(cg_i + 0 + CGROW * 2) += lup1(6);
        tx(cg_i + 1 + CGROW * 2) += lup1(7);
        tx(cg_i + 2 + CGROW * 2) += lup1(8);

        Eigen::Matrix<double, 9, 1> lup2 = scale
            * (DG6_CG2_dX * S12.row(c).transpose() / mesh.hx
                + DG6_CG2_dY * S22.row(c).transpose() / mesh.hy);
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
    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<2>& tx, CGVector<2>& ty, const CellVector<8>& S11,
        const CellVector<8>& S12, const CellVector<8>& S22) const
    {
        const size_t CGROW = 2 * mesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 9, 1> lup1 = scale
            * (DG8_CG2_dX * S11.row(c).transpose() / mesh.hx
                + DG8_CG2_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 2) += lup1(2);
        tx(cg_i + 0 + CGROW) += lup1(3);
        tx(cg_i + 1 + CGROW) += lup1(4);
        tx(cg_i + 2 + CGROW) += lup1(5);
        tx(cg_i + 0 + CGROW * 2) += lup1(6);
        tx(cg_i + 1 + CGROW * 2) += lup1(7);
        tx(cg_i + 2 + CGROW * 2) += lup1(8);

        Eigen::Matrix<double, 9, 1> lup2 = scale
            * (DG8_CG2_dX * S12.row(c).transpose() / mesh.hx
                + DG8_CG2_dY * S22.row(c).transpose() / mesh.hy);
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
    void AddStressTensorCell(const Mesh& mesh, const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<2>& tx, CGVector<2>& ty, const CellVector<1>& S11,
        const CellVector<1>& S12, const CellVector<1>& S22) const
    {
        const size_t CGROW = 2 * mesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        Eigen::Matrix<double, 9, 1> lup1 = scale
            * (DG1_CG2_dX * S11.row(c).transpose() / mesh.hx
                + DG1_CG2_dY * S12.row(c).transpose() / mesh.hy);
        tx(cg_i + 0) += lup1(0);
        tx(cg_i + 1) += lup1(1);
        tx(cg_i + 2) += lup1(2);
        tx(cg_i + 0 + CGROW) += lup1(3);
        tx(cg_i + 1 + CGROW) += lup1(4);
        tx(cg_i + 2 + CGROW) += lup1(5);
        tx(cg_i + 0 + CGROW * 2) += lup1(6);
        tx(cg_i + 1 + CGROW * 2) += lup1(7);
        tx(cg_i + 2 + CGROW * 2) += lup1(8);

        Eigen::Matrix<double, 9, 1> lup2 = scale
            * (DG1_CG2_dX * S12.row(c).transpose() / mesh.hx
                + DG1_CG2_dY * S22.row(c).transpose() / mesh.hy);
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

    void AddStressTensorCell(const SasipMesh& smesh, const double scale, const size_t eid, const size_t cx,
        const size_t cy, CGVector<2>& tmpx, CGVector<2>& tmpy, const CellVector<8>& S11,
        const CellVector<8>& S12, const CellVector<8>& S22) const
    {
        //      (Mv)_i = (v, phi_i) = - (S, nabla Phi_i)

        // (M vx)_i = (vx, phi_i) = - (S11, d_x phi_i) - (S12, d_y phi_i)

        const size_t CGROW = 2 * smesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        const Eigen::Matrix<Nextsim::FloatType, 1, 9> S11_g = S11.row(eid) * PSI<8, 3>; //!< velocity in GP
        const Eigen::Matrix<Nextsim::FloatType, 1, 9> S22_g = S22.row(eid) * PSI<8, 3>; //!< velocity in GP
        const Eigen::Matrix<Nextsim::FloatType, 1, 9> S12_g = S12.row(eid) * PSI<8, 3>; //!< velocity in GP

        // J T^{-T}
        const Eigen::Matrix<Nextsim::FloatType, 2, 9> dxT = (ParametricTools::dxT<3>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<3>.array()).matrix();
        const Eigen::Matrix<Nextsim::FloatType, 2, 9> dyT = (ParametricTools::dyT<3>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<3>.array()).matrix();

        const Eigen::Matrix<Nextsim::FloatType, 9, 9> dx_cg2 = PHIx<2, 3>.array().rowwise() * dyT.row(1).array() - PHIy<2, 3>.array().rowwise() * dxT.row(1).array();

        const Eigen::Matrix<Nextsim::FloatType, 9, 9> dy_cg2 = PHIy<2, 3>.array().rowwise() * dxT.row(0).array() - PHIx<2, 3>.array().rowwise() * dyT.row(0).array();

        const Eigen::Matrix<Nextsim::FloatType, 1, 9> tx = dx_cg2 * S11_g.transpose() + dy_cg2 * S12_g.transpose();
        const Eigen::Matrix<Nextsim::FloatType, 1, 9> ty = dx_cg2 * S12_g.transpose() + dy_cg2 * S22_g.transpose();

        tmpx(cg_i + 0) += -tx(0);
        tmpx(cg_i + 1) += -tx(1);
        tmpx(cg_i + 2) += -tx(2);
        tmpx(cg_i + 0 + CGROW) += -tx(3);
        tmpx(cg_i + 1 + CGROW) += -tx(4);
        tmpx(cg_i + 2 + CGROW) += -tx(5);
        tmpx(cg_i + 0 + CGROW * 2) += -tx(6);
        tmpx(cg_i + 1 + CGROW * 2) += -tx(7);
        tmpx(cg_i + 2 + CGROW * 2) += -tx(8);

        tmpy(cg_i + 0) += -ty(0);
        tmpy(cg_i + 1) += -ty(1);
        tmpy(cg_i + 2) += -ty(2);
        tmpy(cg_i + 0 + CGROW) += -ty(3);
        tmpy(cg_i + 1 + CGROW) += -ty(4);
        tmpy(cg_i + 2 + CGROW) += -ty(5);
        tmpy(cg_i + 0 + CGROW * 2) += -ty(6);
        tmpy(cg_i + 1 + CGROW * 2) += -ty(7);
        tmpy(cg_i + 2 + CGROW * 2) += -ty(8);
    }

    void AddStressTensorCell(const ParametricTransformation<2, 8>& ptrans,
        const SasipMesh& smesh, const double scale, const size_t eid, const size_t cx,
        const size_t cy, CGVector<2>& tmpx, CGVector<2>& tmpy, const CellVector<8>& S11,
        const CellVector<8>& S12, const CellVector<8>& S22) const
    {
        const Eigen::Matrix<Nextsim::FloatType, 9, 1> tx = ptrans.divS1[eid] * S11.row(eid).transpose() + ptrans.divS2[eid] * S12.row(eid).transpose();
        const Eigen::Matrix<Nextsim::FloatType, 9, 1> ty = ptrans.divS1[eid] * S12.row(eid).transpose() + ptrans.divS2[eid] * S22.row(eid).transpose();

        const size_t CGROW = 2 * smesh.nx + 1;
        const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

        tmpx.block<3, 1>(cg_i, 0) -= tx.block<3, 1>(0, 0);
        tmpx.block<3, 1>(cg_i + CGROW, 0) -= tx.block<3, 1>(3, 0);
        tmpx.block<3, 1>(cg_i + 2 * CGROW, 0) -= tx.block<3, 1>(6, 0);

        tmpy.block<3, 1>(cg_i, 0) -= ty.block<3, 1>(0, 0);
        tmpy.block<3, 1>(cg_i + CGROW, 0) -= ty.block<3, 1>(3, 0);
        tmpy.block<3, 1>(cg_i + 2 * CGROW, 0) -= ty.block<3, 1>(6, 0);
    }

    //! Sets the velocity vector to zero along the boundary
    template <int CG>
    void DirichletZero(const Mesh& mesh, CGVector<CG>& v) const;
    //! Sets the velocity vector to zero along the boundary
    template <int CG>
    void DirichletZero(const SasipMesh& smesh, CGVector<CG>& v) const;

    //! Sets the velocity vector for compresion testcase along the boundary
    template <int CG>
    void DirichletCompressionVx(const Mesh& mesh, CGVector<CG>& vector) const;
    template <int CG>
    void DirichletCompressionVy(const Mesh& mesh, CGVector<CG>& vector, const double& value) const;

    template <int DG>
    void DirichletCompressionBottom(const Mesh& mesh, CellVector<DG>& vector, const double& value) const;
    template <int CG>
    void DirichletCompressionBottom(const Mesh& mesh, CGVector<CG>& vector, const double& value) const;

    template <int DG>
    void DirichletCompressionTop(const Mesh& mesh, CellVector<DG>& vector, const double& value) const;
    template <int CG>
    void DirichletCompressionTop(const Mesh& mesh, CGVector<CG>& vector, const double& value) const;

    template <int CG>
    void DirichletCompressionLeft(const Mesh& mesh, CGVector<CG>& vector, const double& value) const;
    template <int CG>
    void DirichletCompressionRight(const Mesh& mesh, CGVector<CG>& vector, const double& value) const;
    template <int CG>
    void DirichletCompressionFixCorner(const Mesh& mesh, CGVector<CG>& vector) const;

    //! Interpolates a DG-Vector to a CG-Vector (OLD MESH)
    template <int CG, int DG>
    void InterpolateDGToCG(const Mesh& mesh, CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    template <int CG, int DG>
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& cg_A, const CellVector<1>& A) const
    {
        const size_t CGDofsPerRow = mesh.nx + 1;
        const size_t cgi
            = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * A(c);
        cg_A(cgi + 1) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow + 1) += 0.25 * A(c);
    }
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& cg_A, const CellVector<3>& A) const
    {
        const size_t CGDofsPerRow = mesh.nx + 1;
        const size_t cgi
            = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow + 1) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    }
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<1>& A) const
    {
        const size_t CGDofsPerRow = 2 * mesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
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
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<3>& A) const
    {
        const size_t CGDofsPerRow = 2 * mesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.5 * (A(c, 0) - 0.5 * A(c, 2));
        cg_A(cgi + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.5 * (A(c, 0) - 0.5 * A(c, 1));
        cg_A(cgi + CGDofsPerRow + 1) += A(c, 0);
        cg_A(cgi + CGDofsPerRow + 2) += 0.5 * (A(c, 0) + 0.5 * A(c, 1));
        cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * (A(c, 0) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    }
    void InterpolateDGToCGCell(const Mesh& mesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<6>& A) const
    {
        const size_t CGDofsPerRow = 2 * mesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.5 * (A(c, 0) - 0.5 * A(c, 2));
        cg_A(cgi + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.5 * (A(c, 0) - 0.5 * A(c, 1));
        cg_A(cgi + CGDofsPerRow + 1) += A(c, 0);
        cg_A(cgi + CGDofsPerRow + 2) += 0.5 * (A(c, 0) + 0.5 * A(c, 1));
        cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * (A(c, 0) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
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

    //! Interpolates a DG-Vector to a CG-Vector (Sasip-Mesh)
    // IS THIS ALL OK and does not depend on the mesh degeneration?
    template <int CG, int DG>
    void InterpolateDGToCG(const SasipMesh& smesh, CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    template <int CG, int DG>
    void InterpolateDGToCGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    void InterpolateDGToCGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& cg_A, const CellVector<1>& A) const
    {
        const size_t CGDofsPerRow = smesh.nx + 1;
        const size_t cgi
            = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * A(c);
        cg_A(cgi + 1) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow) += 0.25 * A(c);
        cg_A(cgi + CGDofsPerRow + 1) += 0.25 * A(c);
    }
    void InterpolateDGToCGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& cg_A, const CellVector<3>& A) const
    {
        const size_t CGDofsPerRow = smesh.nx + 1;
        const size_t cgi
            = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow + 1) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    }
    void InterpolateDGToCGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<1>& A) const
    {
        const size_t CGDofsPerRow = 2 * smesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
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
    void InterpolateDGToCGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<3>& A) const
    {
        const size_t CGDofsPerRow = 2 * smesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.5 * (A(c, 0) - 0.5 * A(c, 2));
        cg_A(cgi + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.5 * (A(c, 0) - 0.5 * A(c, 1));
        cg_A(cgi + CGDofsPerRow + 1) += A(c, 0);
        cg_A(cgi + CGDofsPerRow + 2) += 0.5 * (A(c, 0) + 0.5 * A(c, 1));
        cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * (A(c, 0) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    }
    void InterpolateDGToCGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<6>& A) const
    {
        const size_t CGDofsPerRow = 2 * smesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
        cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + 1) += 0.5 * (A(c, 0) - 0.5 * A(c, 2));
        cg_A(cgi + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
        cg_A(cgi + CGDofsPerRow) += 0.5 * (A(c, 0) - 0.5 * A(c, 1));
        cg_A(cgi + CGDofsPerRow + 1) += A(c, 0);
        cg_A(cgi + CGDofsPerRow + 2) += 0.5 * (A(c, 0) + 0.5 * A(c, 1));
        cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * (A(c, 0) + 0.5 * A(c, 2));
        cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    }

    //! Adjusts the interpolation on the boundary
    template <int CG>
    void InterpolateDGToCGBoundary(const SasipMesh& smesh, CGVector<CG>& cg_A) const;
    void InterpolateDGToCGBoundary(const SasipMesh& smesh, CGVector<1>& cg_A) const
    {
        const size_t CGDofsPerRow = smesh.nx + 1;
        const size_t UpperLeftIndex = CGDofsPerRow * smesh.ny;

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nx + 1; ++i) {
            cg_A(i) *= 2.0;
            cg_A(UpperLeftIndex + i) *= 2.0;
        }
#pragma omp parallel for
        for (size_t i = 0; i < smesh.ny + 1; ++i) {
            cg_A(i * CGDofsPerRow) *= 2.0;
            cg_A(i * CGDofsPerRow + smesh.nx) *= 2.0;
        }
    }
    void InterpolateDGToCGBoundary(const SasipMesh& smesh, CGVector<2>& cg_A) const
    {
        const size_t CGDofsPerRow = 2 * smesh.nx + 1;
        const size_t UpperLeftIndex = 2 * CGDofsPerRow * smesh.ny;

#pragma omp parallel for
        for (size_t i = 0; i < 2 * smesh.nx + 1; ++i) {
            cg_A(i) *= 2.0;
            cg_A(UpperLeftIndex + i) *= 2.0;
        }
#pragma omp parallel for
        for (size_t i = 0; i < 2 * smesh.ny + 1; ++i) {
            cg_A(i * CGDofsPerRow) *= 2.0;
            cg_A(i * CGDofsPerRow + 2 * smesh.nx) *= 2.0;
        }
    }
};
} /* namespace Nextsim */

#endif /* __CGMOMENTUM_HPP */
