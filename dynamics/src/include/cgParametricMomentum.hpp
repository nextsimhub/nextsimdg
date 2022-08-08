/*!
 * @file cgParametricMomentum.hpp
 * @date 5 August 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __CGPARAMETRICMOMENTUM_HPP
#define __CGPARAMETRICMOMENTUM_HPP

#include "ParametricTools.hpp"
#include "cgVector.hpp"
#include "codeGenerationCGToDG.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

template <int CG, int DGstress>
class CGParametricMomentum {
private:
    const SasipMesh& smesh; //!< const-reference to the mesh

    /*!
     * Stores precomputed values for efficient numerics on transformed mesh
     * accelerates numerics but substantial memory effort!
     */
    static constexpr int precompute_matrices = 1;
    ParametricTransformation<CG, DGstress> ptrans;

    //! vectors storing the velocity (node-wise)
    CGVector<CG> vx, vy;

    //! vectors storing ocean and atm velocity (node-wise)
    CGVector<CG> ox, oy, ax, ay;

    //! temporary vectors. Maybe we can remove them?
    CGVector<CG> tmpx, tmpy;

    //! Vector to store the lumpes mass matrix. Is directly initialized when the mesh is known
    CGVector<CG> lumpedcgmass;

    //! Vector to store the CG-Version of ice concentration and ice height
    CGVector<CG> cg_A, cg_H;

public:
    CGParametricMomentum(const SasipMesh& sm)
        : smesh(sm)
    {
        if (!(smesh.nelements > 0)) {
            std::cerr << "CGParametricMomentum: The mesh has to be initialized first!" << std::endl;
            abort();
        }

        // just simple sanity checks
        assert(smesh.nelements == smesh.nx * smesh.ny);
        assert(smesh.nnodes == (smesh.nx + 1) * (smesh.ny + 1));

        // resize the vectors
        vx.resize_by_mesh(smesh);
        vy.resize_by_mesh(smesh);
        vx.setZero();
        vy.setZero();

        cg_A.resize_by_mesh(smesh);
        cg_H.resize_by_mesh(smesh);

        ax.resize_by_mesh(smesh);
        ay.resize_by_mesh(smesh);
        ox.resize_by_mesh(smesh);
        oy.resize_by_mesh(smesh);

        tmpx.resize_by_mesh(smesh);
        tmpy.resize_by_mesh(smesh);

        //! precomputes matrices
        if (precompute_matrices == 1)
            ptrans.BasicInit(smesh);

        // initialize the lumped mass
        lumpedcgmass.resize_by_mesh(smesh);
        ParametricTools::lumpedCGMassMatrix(smesh, lumpedcgmass);
    }

    // Access to members
    const CGVector<CG>& GetVx() const { return vx; }
    const CGVector<CG>& GetVy() const { return vy; }
    CGVector<CG>& GetVx() { return vx; }
    CGVector<CG>& GetVy() { return vy; }

    const CGVector<CG>& GetOceanx() const { return ox; }
    const CGVector<CG>& GetOceany() const { return oy; }
    CGVector<CG>& GetOceanx() { return ox; }
    CGVector<CG>& GetOceany() { return oy; }
    const CGVector<CG>& GetAtmx() const { return ax; }
    const CGVector<CG>& GetAtmy() const { return ay; }
    CGVector<CG>& GetAtmx() { return ax; }
    CGVector<CG>& GetAtmy() { return ay; }

    // High level Functions

    //! performs one complete mEVP cycle with NT_evp subiterations
    template <int DG>
    void mEVPIteration(size_t NT_evp, double alpha, double beta,
        double dt_adv,
        const CellVector<DG>& H, const CellVector<DG>& A,
        CellVector<DGstress>& E11, CellVector<DGstress>& E12, CellVector<DGstress>& E22,
        CellVector<DGstress>& S11, CellVector<DGstress>& S12, CellVector<DGstress>& S22);

    /*!
     * The following functions take care of the interpolation and projection
     * between CG and DG functions
     */
    //! Projects a CG function to a DG function
    template <int DG>
    void ProjectCGToDG(CellVector<DG>& dg, const CGVector<CG>& cg);

    //! Projects the symmetric gradient of the CG2 velocity into the DG1 space
    template <int DG>
    void ProjectCG2VelocityToDG1Strain(CellVector<DG>& E11, CellVector<DG>& E12,
        CellVector<DG>& E22);

    /*!
     * Evaluates (S, nabla phi) and adds it to tx/ty - Vector
     */
    template <int DG>
    void AddStressTensor(const double scale, CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    template <int DG>
    void AddStressTensorCell(const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<CG>& tx, CGVector<CG>& ty, const CellVector<DG>& S11,
        const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    void AddStressTensorCell(const double scale, const size_t eid, const size_t cx,
        const size_t cy, CGVector<2>& tmpx, CGVector<2>& tmpy, const CellVector<8>& S11,
        const CellVector<8>& S12, const CellVector<8>& S22) const
    {
        if (precompute_matrices == 0) // all is compute on-the-fly
        {
            //      (Mv)_i = (v, phi_i) = - (S, nabla Phi_i)

            // (M vx)_i = (vx, phi_i) = - (S11, d_x phi_i) - (S12, d_y phi_i)

            const size_t CGROW = 2 * smesh.nx + 1;
            const size_t cg_i = 2 * CGROW * cy + 2 * cx; //!< lower left CG-index in element (cx,cy)

            const Eigen::Matrix<Nextsim::FloatType, 1, 9> S11_g = S11.row(eid) * BiG83; //!< velocity in GP
            const Eigen::Matrix<Nextsim::FloatType, 1, 9> S22_g = S22.row(eid) * BiG83; //!< velocity in GP
            const Eigen::Matrix<Nextsim::FloatType, 1, 9> S12_g = S12.row(eid) * BiG83; //!< velocity in GP

            // J T^{-T}
            const Eigen::Matrix<Nextsim::FloatType, 2, 9> dxT = (ParametricTools::dxT<3>(smesh, eid).array().rowwise() * GAUSSWEIGHTS_3.array()).matrix();
            const Eigen::Matrix<Nextsim::FloatType, 2, 9> dyT = (ParametricTools::dyT<3>(smesh, eid).array().rowwise() * GAUSSWEIGHTS_3.array()).matrix();

            const Eigen::Matrix<Nextsim::FloatType, 9, 9> dx_cg2 = CG_CG2_dx_in_GAUSS3.array().rowwise() * dyT.row(1).array() - CG_CG2_dy_in_GAUSS3.array().rowwise() * dxT.row(1).array();

            const Eigen::Matrix<Nextsim::FloatType, 9, 9> dy_cg2 = CG_CG2_dy_in_GAUSS3.array().rowwise() * dxT.row(0).array() - CG_CG2_dx_in_GAUSS3.array().rowwise() * dyT.row(0).array();

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
        } else if (precompute_matrices == 1) // use precomputed values
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
        } else
            abort();
    }

    //! Sets the velocity vector to zero along the boundary
    void DirichletZero()
    {
        DirichletZero(vx);
        DirichletZero(vy);
    }
    void DirichletZero(CGVector<CG>& v);

    //! Interpolates a DG-Vector to a CG-Vector (Sasip-Mesh)
    // IS THIS ALL OK and does not depend on the mesh degeneration?
    template <int DG>
    void InterpolateDGToCG(CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    template <int DG>
    void InterpolateDGToCGCell(const size_t c, const size_t cx, const size_t cy,
        CGVector<CG>& cg_A, const CellVector<DG>& A) const;

    void InterpolateDGToCGCell(const size_t c, const size_t cx, const size_t cy,
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
    void InterpolateDGToCGCell(const size_t c, const size_t cx, const size_t cy,
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
    void InterpolateDGToCGCell(const size_t c, const size_t cx, const size_t cy,
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
    void InterpolateDGToCGCell(const size_t c, const size_t cx, const size_t cy,
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
    void InterpolateDGToCGCell(const size_t c, const size_t cx, const size_t cy,
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
    void InterpolateDGToCGBoundary(CGVector<CG>& cg_A) const;
};

template <>
void CGParametricMomentum<1, 3>::InterpolateDGToCGBoundary(CGVector<1>& cg_A) const
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

template <>
void CGParametricMomentum<2, 8>::InterpolateDGToCGBoundary(CGVector<2>& cg_A) const
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

} /* namespace Nextsim */

#endif /* __CGMOMENTUM_HPP */
