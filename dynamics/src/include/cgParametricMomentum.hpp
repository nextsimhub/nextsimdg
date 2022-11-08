/*!
 * @file cgParametricMomentum.hpp
 * @date 5 August 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __CGPARAMETRICMOMENTUM_HPP
#define __CGPARAMETRICMOMENTUM_HPP

#include "MEBParameters.hpp"
#include "ParametricTools.hpp"
#include "VPParameters.hpp"
#include "cgVector.hpp"
#include "codeGenerationCGToDG.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

template <int CG, int DGstress>
class CGParametricMomentum {
private:
    const ParametricMesh& smesh; //!< const-reference to the mesh

    /*!
     * Stores precomputed values for efficient numerics on transformed mesh
     * accelerates numerics but substantial memory effort!
     */
    static constexpr int precompute_matrices = 1;
    ParametricTransformation<CG, DGstress> ptrans;

    //! vectors storing the velocity (node-wise)
    CGVector<CG> vx, vy;

    //! vectors storing the average velocity in sub-iterations (node-wise)
    CGVector<CG> avg_vx, avg_vy;

    //! vectors storing ocean and atm velocity (node-wise)
    CGVector<CG> ox, oy, ax, ay;

    ////// Internal variables

    //! temporary vectors. Maybe we can remove them?
    CGVector<CG> tmpx, tmpy;

    //! old velocities. They are required during temporary vectors. Maybe we can remove them?
    CGVector<CG> vx_mevp, vy_mevp;

    //! Vector to store the lumpes mass matrix. Is directly initialized when the mesh is known
    CGVector<CG> lumpedcgmass;

    //! Vector to store the CG-Version of ice concentration and ice height
    CGVector<CG> cg_A, cg_H, cg_D;

public:
    //! Vectors storing strain and strss
    DGVector<DGstress> E11, E12, E22;
    DGVector<DGstress> S11, S12, S22;

public:
    CGParametricMomentum(const ParametricMesh& sm)
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

        avg_vx.resize_by_mesh(smesh);
        avg_vy.resize_by_mesh(smesh);
        avg_vx.setZero();
        avg_vy.setZero();

        cg_A.resize_by_mesh(smesh);
        cg_H.resize_by_mesh(smesh);
        cg_D.resize_by_mesh(smesh);

        ax.resize_by_mesh(smesh);
        ay.resize_by_mesh(smesh);
        ox.resize_by_mesh(smesh);
        oy.resize_by_mesh(smesh);

        tmpx.resize_by_mesh(smesh);
        tmpy.resize_by_mesh(smesh);

        E11.resize_by_mesh(smesh);
        E12.resize_by_mesh(smesh);
        E22.resize_by_mesh(smesh);
        S11.resize_by_mesh(smesh);
        S12.resize_by_mesh(smesh);
        S22.resize_by_mesh(smesh);

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

    const CGVector<CG>& GetAvgSubiterVx() const { return avg_vx; }
    const CGVector<CG>& GetAvgSubiterVy() const { return avg_vy; }
    CGVector<CG>& GetAvgSubiterVx() { return avg_vx; }
    CGVector<CG>& GetAvgSubiterVy() { return avg_vy; }

    const CGVector<CG>& GetOceanx() const { return ox; }
    const CGVector<CG>& GetOceany() const { return oy; }
    CGVector<CG>& GetOceanx() { return ox; }
    CGVector<CG>& GetOceany() { return oy; }
    const CGVector<CG>& GetAtmx() const { return ax; }
    const CGVector<CG>& GetAtmy() const { return ay; }
    CGVector<CG>& GetAtmx() { return ax; }
    CGVector<CG>& GetAtmy() { return ay; }

    const DGVector<DGstress> GetE11() const { return E11; }
    const DGVector<DGstress> GetE12() const { return E12; }
    const DGVector<DGstress> GetE22() const { return E22; }

    const DGVector<DGstress> GetS11() const { return S11; }
    const DGVector<DGstress> GetS12() const { return S12; }
    const DGVector<DGstress> GetS22() const { return S22; }

    // High level Functions

    /*!
     *  prepare the subcucling iteration:
     *  - store old velocity
     *  - interpoalte ice height & concentration ( & damage) to cg
     */
    template <int DG>
    void prepareIteration(const DGVector<DG>& H, const DGVector<DG>& A);
    template <int DG>
    void prepareIteration(const DGVector<DG>& H, const DGVector<DG>& A,
        const DGVector<DG>& D);

    //! performs one complete mEVP cycle with NT_evp subiterations
    template <int DG>
    void mEVPStep(const VPParameters& vpparameters,
        size_t NT_evp, double alpha, double beta,
        double dt_adv,
        const DGVector<DG>& H, const DGVector<DG>& A);

    //! performs one complete MEB timestep with NT_meb subiterations
    template <int DG>
    void MEBStep(const MEBParameters& vpparameters, size_t NT_meb,
        double dt_adv, const DGVector<DG>& H, const DGVector<DG>& A, DGVector<DG>& D);

    //! performs one complete MEB timestep with NT_meb subiterations
    template <int DG>
    void MEBIteration(const MEBParameters& vpparameters, size_t NT_meb,
        double dt_adv, const DGVector<DG>& H, const DGVector<DG>& A, DGVector<DG>& D);

    /*!
     * The following functions take care of the interpolation and projection
     * between CG and DG functions
     */
    //! Projects the symmetric gradient of the CG velocity into the DG space
    void ProjectCGVelocityToDGStrain();

    /*!
     * Evaluates (S, nabla phi) and adds it to tx/ty - Vector
     */
    void AddStressTensor(const double scale, CGVector<CG>& tx, CGVector<CG>& ty) const;

    void AddStressTensorCell(const double scale, const size_t c, const size_t cx,
        const size_t cy, CGVector<CG>& tx, CGVector<CG>& ty) const;

    //! Sets the velocity vector to zero along the boundary
    void DirichletZero()
    {
        DirichletZero(vx);
        DirichletZero(vy);
    }
    void DirichletZero(CGVector<CG>& v);
};

template <int CG, int DG>
void CGParametricMomentum<CG, DG>::AddStressTensorCell(const double scale, const size_t eid, const size_t cx,
    const size_t cy, CGVector<CG>& tmpx, CGVector<CG>& tmpy) const
{
#define NGP (DG == 8 ? 3 : (DG == 3 ? 2 : 1))
    if (precompute_matrices == 0) // all is compute on-the-fly
    {
        //      (Mv)_i = (v, phi_i) = - (S, nabla Phi_i)

        // (M vx)_i = (vx, phi_i) = - (S11, d_x phi_i) - (S12, d_y phi_i)

        const size_t CGROW = CG * smesh.nx + 1;
        const size_t cg_i = CG * CGROW * cy + CG * cx; //!< lower left CG-index in element (cx,cy)

        const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> S11_g = S11.row(eid) * PSI<DG, NGP>; //!< stress in GP
        const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> S22_g = S22.row(eid) * PSI<DG, NGP>;
        const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> S12_g = S12.row(eid) * PSI<DG, NGP>;

        // J T^{-T}
        const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dxT = (ParametricTools::dxT<NGP>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<NGP>.array()).matrix();
        const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dyT = (ParametricTools::dyT<NGP>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<NGP>.array()).matrix();

        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), NGP* NGP> dx_cg2 = PHIx<CG, NGP>.array().rowwise() * dyT.row(1).array() - PHIy<CG, NGP>.array().rowwise() * dxT.row(1).array();

        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), NGP* NGP> dy_cg2 = PHIy<CG, NGP>.array().rowwise() * dxT.row(0).array() - PHIx<CG, NGP>.array().rowwise() * dyT.row(0).array();

        const Eigen::Matrix<Nextsim::FloatType, 1, (CG == 2 ? 9 : 4)> tx = dx_cg2 * S11_g.transpose() + dy_cg2 * S12_g.transpose();
        const Eigen::Matrix<Nextsim::FloatType, 1, (CG == 2 ? 9 : 4)> ty = dx_cg2 * S12_g.transpose() + dy_cg2 * S22_g.transpose();

        if (CG == 1) {
            tmpx(cg_i + 0) += -tx(0);
            tmpx(cg_i + 1) += -tx(1);
            tmpx(cg_i + 0 + CGROW) += -tx(2);
            tmpx(cg_i + 1 + CGROW) += -tx(3);

            tmpy(cg_i + 0) += -ty(0);
            tmpy(cg_i + 1) += -ty(1);
            tmpy(cg_i + 0 + CGROW) += -ty(2);
            tmpy(cg_i + 1 + CGROW) += -ty(3);
        } else if (CG == 2) {
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
    } else if (precompute_matrices == 1) // use precomputed values
    {
        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), 1> tx = ptrans.divS1[eid] * S11.row(eid).transpose() + ptrans.divS2[eid] * S12.row(eid).transpose();
        const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), 1> ty = ptrans.divS1[eid] * S12.row(eid).transpose() + ptrans.divS2[eid] * S22.row(eid).transpose();

        const size_t CGROW = CG * smesh.nx + 1;
        const size_t cg_i = CG * CGROW * cy + CG * cx; //!< lower left CG-index in element (cx,cy)

        if (CG == 1) {
            tmpx(cg_i + 0) += -tx(0);
            tmpx(cg_i + 1) += -tx(1);
            tmpx(cg_i + 0 + CGROW) += -tx(2);
            tmpx(cg_i + 1 + CGROW) += -tx(3);

            tmpy(cg_i + 0) += -ty(0);
            tmpy(cg_i + 1) += -ty(1);
            tmpy(cg_i + 0 + CGROW) += -ty(2);
            tmpy(cg_i + 1 + CGROW) += -ty(3);
        } else if (CG == 2) {
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

            // tmpx.block<CG+1, 1>(cg_i, 0) -= tx.block<CG+1, 1>(0, 0);
            // tmpx.block<CG+1, 1>(cg_i + CGROW, 0) -= tx.block<CG+1, 1>(CG+1, 0);
            // tmpx.block<CG+1, 1>(cg_i + 2 * CGROW, 0) -= tx.block<CG+1, 1>(2*CG+2, 0);

            // tmpy.block<CG+1, 1>(cg_i, 0) -= ty.block<CG+1, 1>(0, 0);
            // tmpy.block<CG+1, 1>(cg_i + CGROW, 0) -= ty.block<CG+1, 1>(CG+1, 0);
            // tmpy.block<CG+1, 1>(cg_i + 2 * CGROW, 0) -= ty.block<CG+1, 1>(2*CG+2, 0);
        }
    } else
        abort();

#undef NGP
}

} /* namespace Nextsim */

#endif /* __CGMOMENTUM_HPP */
