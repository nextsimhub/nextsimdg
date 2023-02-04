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
#include "ParametricMap.hpp"

namespace Nextsim {

  
template <int CG>
class CGParametricMomentum {
private:
    const ParametricMesh& smesh; //!< const-reference to the mesh

  const COORDINATES coordinatesystem; //! const-ref to the coordinate system (2d Cartesian / SPHERICAL / ...)
  
    /*!
     * Stores precomputed values for efficient numerics on transformed mesh
     * accelerates numerics but substantial memory effort!
     */
    static constexpr int precompute_matrices = 1;

  ParametricMomentumMap<CG> pmap;

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


    //! Vector to store the CG-Version of ice concentration and ice height
    CGVector<CG> cg_A, cg_H, cg_D;

    //! Vectors storing strain and strss
    DGVector<CG2DGSTRESS(CG)> E11, E12, E22;
    DGVector<CG2DGSTRESS(CG)> S11, S12, S22;

public:
  CGParametricMomentum(const ParametricMesh& sm, const COORDINATES coords)
    : smesh(sm), coordinatesystem(coords), pmap(sm,coords)
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


        /*!
	 * initialize the lumped mass
	 * At periodic boundaries, the values must be added from both sides
	 */
	pmap.InitializeLumpedCGMassMatrix();

	/*!
	 * Compute matrices for performing mEVP / BBM / MEB etc. stress updates
	 */
	pmap.InitializeDivSMatrices();
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

    const DGVector<CG2DGSTRESS(CG)> GetE11() const { return E11; }
    const DGVector<CG2DGSTRESS(CG)> GetE12() const { return E12; }
    const DGVector<CG2DGSTRESS(CG)> GetE22() const { return E22; }

    const DGVector<CG2DGSTRESS(CG)> GetS11() const { return S11; }
    const DGVector<CG2DGSTRESS(CG)> GetS12() const { return S12; }
    const DGVector<CG2DGSTRESS(CG)> GetS22() const { return S22; }
  
  const CGVector<CG>& GetcgH() const { return cg_H; }
  const CGVector<CG>& GetcgA() const { return cg_A; }
  const CGVector<CG>& GetcgD() const { return cg_D; }
   CGVector<CG>& GetcgH()  { return cg_H; }
   CGVector<CG>& GetcgA()  { return cg_A; }
   CGVector<CG>& GetcgD()  { return cg_D; }
  
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

  /*!
   * AddPeriodic is to be called, after (sigma, Nabla Phi) is computed
   * On periodic boundaries, the contributions from both sides must be added
   */
  void AddPeriodic(CGVector<CG>& v);
  /*!
   * AveragePeriodic replaces the values on both sides by
   * the average of them
   */
  void AveragePeriodic(CGVector<CG>& v);
  void CheckPeriodicity(CGVector<CG>& v);
};

template <int CG>
void CGParametricMomentum<CG>::AddStressTensorCell(const double scale, const size_t eid, const size_t cx,
    const size_t cy, CGVector<CG>& tmpx, CGVector<CG>& tmpy) const
{
  // pick the number of Gauss points according to the degree
  
#define NGP (CG == 1 ? 2 : 3) 
    // if (precompute_matrices == 0) // all is compute on-the-fly
    // {
    //     //      (Mv)_i = (v, phi_i) = - (S, nabla Phi_i)

    //     // (M vx)_i = (vx, phi_i) = - (S11, d_x phi_i) - (S12, d_y phi_i)

    //     const size_t CGROW = CG * smesh.nx + 1;
    //     const size_t cg_i = CG * CGROW * cy + CG * cx; //!< lower left CG-index in element (cx,cy)

    //     const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> S11_g = S11.row(eid) * PSI<CG2DGSTRESS(CG), NGP>; //!< stress in GP
    //     const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> S22_g = S22.row(eid) * PSI<CG2DGSTRESS(CG), NGP>;
    //     const Eigen::Matrix<Nextsim::FloatType, 1, NGP* NGP> S12_g = S12.row(eid) * PSI<CG2DGSTRESS(CG), NGP>;

    //     // J T^{-T}
    //     const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dxT = (ParametricTools::dxT<NGP>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<NGP>.array()).matrix();
    //     const Eigen::Matrix<Nextsim::FloatType, 2, NGP* NGP> dyT = (ParametricTools::dyT<NGP>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<NGP>.array()).matrix();

    //     const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), NGP* NGP> dx_cg2 = PHIx<CG, NGP>.array().rowwise() * dyT.row(1).array() - PHIy<CG, NGP>.array().rowwise() * dxT.row(1).array();

    //     const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), NGP* NGP> dy_cg2 = PHIy<CG, NGP>.array().rowwise() * dxT.row(0).array() - PHIx<CG, NGP>.array().rowwise() * dyT.row(0).array();

    //     const Eigen::Matrix<Nextsim::FloatType, 1, (CG == 2 ? 9 : 4)> tx = scale * (dx_cg2 * S11_g.transpose() + dy_cg2 * S12_g.transpose());
    //     const Eigen::Matrix<Nextsim::FloatType, 1, (CG == 2 ? 9 : 4)> ty = scale * (dx_cg2 * S12_g.transpose() + dy_cg2 * S22_g.transpose());

    //     if (CG == 1) {
    //         tmpx(cg_i + 0) += -tx(0);
    //         tmpx(cg_i + 1) += -tx(1);
    //         tmpx(cg_i + 0 + CGROW) += -tx(2);
    //         tmpx(cg_i + 1 + CGROW) += -tx(3);

    //         tmpy(cg_i + 0) += -ty(0);
    //         tmpy(cg_i + 1) += -ty(1);
    //         tmpy(cg_i + 0 + CGROW) += -ty(2);
    //         tmpy(cg_i + 1 + CGROW) += -ty(3);
    //     } else if (CG == 2) {
    //         tmpx(cg_i + 0) += -tx(0);
    //         tmpx(cg_i + 1) += -tx(1);
    //         tmpx(cg_i + 2) += -tx(2);
    //         tmpx(cg_i + 0 + CGROW) += -tx(3);
    //         tmpx(cg_i + 1 + CGROW) += -tx(4);
    //         tmpx(cg_i + 2 + CGROW) += -tx(5);
    //         tmpx(cg_i + 0 + CGROW * 2) += -tx(6);
    //         tmpx(cg_i + 1 + CGROW * 2) += -tx(7);
    //         tmpx(cg_i + 2 + CGROW * 2) += -tx(8);

    //         tmpy(cg_i + 0) += -ty(0);
    //         tmpy(cg_i + 1) += -ty(1);
    //         tmpy(cg_i + 2) += -ty(2);
    //         tmpy(cg_i + 0 + CGROW) += -ty(3);
    //         tmpy(cg_i + 1 + CGROW) += -ty(4);
    //         tmpy(cg_i + 2 + CGROW) += -ty(5);
    //         tmpy(cg_i + 0 + CGROW * 2) += -ty(6);
    //         tmpy(cg_i + 1 + CGROW * 2) += -ty(7);
    //         tmpy(cg_i + 2 + CGROW * 2) += -ty(8);
    //     }
    // } else if (precompute_matrices == 1) // use precomputed values
    // {



  
  Eigen::Vector<Nextsim::FloatType, CGDOFS(CG)> tx = scale * (pmap.divS1[eid] * S11.row(eid).transpose() + pmap.divS2[eid] * S12.row(eid).transpose());
  Eigen::Vector<Nextsim::FloatType, CGDOFS(CG)> ty = scale * (pmap.divS1[eid] * S12.row(eid).transpose() + pmap.divS2[eid] * S22.row(eid).transpose());

  if (coordinatesystem == SPHERICAL) // In spherical coordinates there is the additional 'derivative term' arising from the derivative of the units
    {
      ty += scale * pmap.divM[eid] * S11.row(eid).transpose();
      tx += scale * pmap.divM[eid] * S12.row(eid).transpose();
    }
  
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
  }
  else
    abort();


  //   // } else
  //   //     abort();

#undef NGP
}

} /* namespace Nextsim */

#endif /* __CGMOMENTUM_HPP */
