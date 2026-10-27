/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICMAP_HPP
#define __PARAMETRICMAP_HPP

#include "NextsimDynamics.hpp"
#include "ParametricMesh.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"

#include <array>
#include <vector>

namespace Nextsim {

/*!
 * Stores precomputed matrices and stencils that are required
 * for the advection.
 *
 * The coordiante system is encoded into the matrices
 */
template <int DG> class ParametricTransportMap {
    //! Reference to the map. Given with constructor
    const ParametricMesh& smesh;

public:
    //! These terms are required for the cell-term in the advection  -(vA, nabla PHI)
    std::vector<Eigen::Matrix<FloatType, DG, GAUSSPOINTS(DG)>> AdvectionCellTermX;
    std::vector<Eigen::Matrix<FloatType, DG, GAUSSPOINTS(DG)>> AdvectionCellTermY;

    //! The inverse of the dG mass matrix
    std::vector<Eigen::Matrix<FloatType, DG, DG>> InverseDGMassMatrix;

    ParametricTransportMap(const ParametricMesh& sm)
        : smesh(sm)
    {
    }

    //! initialization of the different forms

    //! For cell_term in Advection. Fills AdvectionCellTermX/Y
    void InitializeAdvectionCellTerms();

    /*!
     *
     *  For cell_term in Advection. Fills AdvectionCellTermX/Y
     *
     * In Spherical coordinates the inverse mass is scaled by 1/Radius to account for
     * the proper scaling.
     *
     * R^2 (v',phi) - R (vA, nabla phi) + R <va * N, phi> = 0
     */
    void InitializeInverseDGMassMatrix();
};

/*!
 * Stores precomputed matrices and stencils that are required
 * for the momentum equation
 *
 * The coordiante system is encoded into the matrices
 */
template <int CG, int DG> class ParametricMomentumMap {
    //! Reference to the map. Given with constructor
    const ParametricMesh& smesh;

public:
    //! Vector to store the lumpes mass matrix. Is directly initialized when the mesh is known
    CGVector<CG> lumpedcgmass;
    //! Vector to store the lumpes mass matrix in CG1. Neede to compute SeasurfaceGradient
    CGVector<1> lumpedcg1mass;
    //! Vector to store the landmask in the CG-degrees of freedom
    CGVector<CG> cglandmask;

    /*!
     * These matrices realize the integration of (-div S, phi) = (S, nabla phi)
     * as matrix-vector producs (divS1 * S11 + divS2 * S12 ; divS1 * S21 + divS2 * S22)
     * [ where S12= S21 ]
     * divM is in addition required for spherical coordinates
     */
    using DivMatrix = Eigen::Matrix<FloatType, CGDOFS(CG), CG2DGSTRESS(CG)>;
    std::vector<DivMatrix, Eigen::aligned_allocator<DivMatrix>> divS1, divS2, divM;

    /*!
     * These matrices are used to compute the gradient of the sea surface height via
     * ( gH, Phi) = ( d_[X/Y] SSH, Phi)
     * where SSH is CG1-representation of seaSurfaceHeight
     *
     * Very similar to divS1 and divS2 but working in CG(1) vectors
     */
    using DSSHMatrix = Eigen::Matrix<FloatType, 4, 4>;
    std::vector<DSSHMatrix, Eigen::aligned_allocator<DSSHMatrix>> dX_SSH, dY_SSH;

    /*!
     * These matrices realize the integration of (E, \grad phi) scaled with the
     * inverse mass matrix;
     */
    using GradMatrix = Eigen::Matrix<FloatType, CG2DGSTRESS(CG), CGDOFS(CG)>;
    std::vector<GradMatrix, Eigen::aligned_allocator<GradMatrix>> iMgradX, iMgradY, iMM;

    /*!
     * These matrices are M^-1 J w PSI_i(q),
     * inverse of mass matrix (in DGstress-degree) multiplied with map J, weights w and DG Test
     * functions, all in the GAUSS points, i.e. iMJwPSI \in R^(DGstress x Ngauss)
     */
    using GaussMapMatrix = Eigen::Matrix<FloatType, CG2DGSTRESS(CG), GAUSSPOINTS(CG2DGSTRESS(CG))>;
    std::vector<GaussMapMatrix, Eigen::aligned_allocator<GaussMapMatrix>> iMJwPSI;

    /*!
     * These matrices are M^-1 J w PSI_i(q),
     * but, instead of DG-stress based on the DGAdvection (used for damage)
     */
    using GaussMapAdvectMatrix = Eigen::Matrix<FloatType, DG, GAUSSPOINTS(CG2DGSTRESS(CG))>;
    std::vector<GaussMapAdvectMatrix, Eigen::aligned_allocator<GaussMapAdvectMatrix>> iMJwPSI_dam;

    ParametricMomentumMap(const ParametricMesh& sm)
        : smesh(sm)
    {
    }

    //! initialization of the different forms

    //! initializes the lumped mass matrix for the stress update
    void InitializeLumpedCGMassMatrix();
    //! initializes div-matrices for the stress update
    void InitializeDivSMatrices();
    //! initializes the CG landmask
    void InitializeCGLandmask();
};

/*!
 * Stores precomputed matrices and stencils that are required for the momentum
 * equation with a DG velocity space and the LDG discretization of the stress
 * divergence (DG velocity, DGstress stress, DGstrain strain).
 *
 * Only CARTESIAN coordinates are supported so far.
 */
template <int DG, int DGstress, int DGstrain> class ParametricDGMomentumMap {
    //! Reference to the map. Given with constructor
    const ParametricMesh& smesh;

public:
    template <int R, int C> using Matrix = Eigen::Matrix<FloatType, R, C>;
    template <int R, int C>
    using MatrixVector = std::vector<Matrix<R, C>, Eigen::aligned_allocator<Matrix<R, C>>>;

    /*!
     * The dG mass matrix weighted with a DG field (e.g. rho_ice * H for the inertia term)
     * and its inverse. Depends on the field, InitializeWeightedDGMassMatrix() must be
     * called whenever the field changes.
     */
    MatrixVector<DG, DG> WeightedDGMassMatrix;
    MatrixVector<DG, DG> InverseWeightedDGMassMatrix;

    /*!
     * Cell term (sigma, nabla v) of the velocity equation, used as
     * (dgDivS1 * S11 + dgDivS2 * S12 ; dgDivS1 * S12 + dgDivS2 * S22).
     * No inverse mass included: the velocity system is solved as a whole.
     */
    MatrixVector<DG, DGstress> dgDivS1, dgDivS2;

    /*!
     * Cell term (u, div tau) of the strain equation, multiplied with the local inverse
     * strain mass matrix.
     */
    MatrixVector<DGstrain, DG> iMgradTauX, iMgradTauY;

    //! Local inverse strain mass matrix, also applied to the flux terms of the strain equation
    MatrixVector<DGstrain, DGstrain> InverseStrainMassMatrix;

    /*!
     * M_stress^-1 J w PSI: projection of a field given in the Gauss points onto the
     * DGstress coefficients, as ParametricMomentumMap::iMJwPSI. Used in the stress update.
     */
    MatrixVector<DGstress, GAUSSPOINTS(DGstress)> iMJwPSI;

    /*!
     * J w PSI: projection of a forcing given in the Gauss points onto the (not inverted)
     * load vector of the DG velocity.
     */
    MatrixVector<DG, GAUSSPOINTS(DG)> JwPSI;

    /*!
     * LDG jump penalty J, the part of -<sigma_hat.n, [[v]]> that depends on u, treated
     * implicitly in the velocity system ((1+beta) rho M_H + dt J) u = F:
     *
     *   J(u,v) = sum_{interior e} w_e int_ref [[u]][[v]] dxi
     *          + sum_{Dirichlet e} w_e int_ref u v dxi
     *
     * with the edge weights (edge length l_e, cell size h = sqrt(area))
     *   interior  : w_e = b * l_e / (0.5*(h(c1)+h(c2)))
     *   Dirichlet : w_e = b * l_e / (h(c) * (0.5-a))
     * Side numbering s: 0 bottom, 1 right, 2 top, 3 left; the same edge seen from the
     * neighbour across side s is its side (s+2)%4.
     */

    //! TraceMassSameCell[s](i,j) = int_0^1 phi_i|_s phi_j|_s dxi
    std::array<Matrix<DG, DG>, 4> TraceMassSameCell;

    //! TraceMassNeighborCell[s](i,j) = int_0^1 phi_i|_s phi_j|_(s+2)%4 dxi (trial
    //! function from the neighbour across side s)
    std::array<Matrix<DG, DG>, 4> TraceMassNeighborCell;

    //! Interior edge weight w_e, one value per edge (0 on boundary edges and on edges
    //! touching land). Numbering as ParametricMesh::coordinatesOfEdgeX/Y:
    //!   X-edges, nx*(ny+1): cell c has bottom c, top c+nx
    //!   Y-edges, (nx+1)*ny: cell c=(ix,iy) has left iy*(nx+1)+ix, right +1
    std::vector<FloatType> JumpPenaltyWeightX, JumpPenaltyWeightY;

    //! Diagonal block of J for each cell: sum over its interior and Dirichlet edges of
    //! w_e * TraceMassSameCell[s]
    MatrixVector<DG, DG> JumpPenaltyDiagonalBlock;

    ParametricDGMomentumMap(const ParametricMesh& sm)
        : smesh(sm)
    {
    }

    //! Fills WeightedDGMassMatrix and InverseWeightedDGMassMatrix
    template <int DGweight> void InitializeWeightedDGMassMatrix(const DGVector<DGweight>& weight);

    //! Fills dgDivS1/dgDivS2, iMgradTauX/iMgradTauY and InverseStrainMassMatrix
    void InitializeDGDivMatrices();

    //! Fills JwPSI
    void InitializeJwPSI();

    //! Fills iMJwPSI. Must be called after InitializeDGDivMatrices()
    void InitializeIMJwPSI();

    //! Fills TraceMassSameCell/NeighborCell, JumpPenaltyWeightX/Y and
    //! JumpPenaltyDiagonalBlock for the LDG parameters a, b
    void InitializeJumpPenalty(const FloatType a, const FloatType b);

    /*!
     * y = massCoefficient * M_H x + penaltyCoefficient * J x
     * Uses WeightedDGMassMatrix unless massCoefficient == 0.
     * Velocity operator of the mEVP step: massCoefficient = (1+beta)*rho_ice,
     * penaltyCoefficient = dt.
     */
    void ApplyMassPlusJumpPenalty(const FloatType massCoefficient,
        const FloatType penaltyCoefficient, const DGVector<DG>& x, DGVector<DG>& y) const;
};

}

#endif // #define  __PARAMETRICMAP_HPP
