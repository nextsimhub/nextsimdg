/*!
 * @file LDGMomentum.hpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 * Local discontinuous Galerkin (LDG) discretization of the momentum equation with the mEVP
 * solver. Velocity, strain and stress are in the same DG space (DG components), ice height
 * and concentration in the advection space (DGadvection components).
 */

#ifndef __LDGMOMENTUM_HPP
#define __LDGMOMENTUM_HPP

#include "ParametricMap.hpp"
#include "ParametricMesh.hpp"
#include "VPParameters.hpp"
#include "dgVector.hpp"

#include <array>
#include <vector>

namespace Nextsim {

/*!
 * LDG edge term -<sigma_hat.n, [[v]]> of the velocity equation for a given stress, ADDED to
 * vx_residual/vy_residual (not mass-inverted).
 *   interior edges  : sigma_hat.n = (0.5-a) sigma1.n + (0.5+a) sigma2.n
 *   Dirichlet edges : sigma_hat.n = sigma.n (smesh.dirichlet)
 *   free edges      : sigma_hat.n = 0 (outer mesh boundary without land), no contribution
 * The jump penalty -b[[u]] of sigma_hat.n is not included, it acts on the unknown u and is
 * part of the velocity operator (ParametricDGMomentumMap::ApplyMassPlusJumpPenalty).
 */
template <int DG, int DGstress>
void ldgSigmaFlux(const ParametricMesh& smesh, const FloatType a, const DGVector<DGstress>& S11,
    const DGVector<DGstress>& S12, const DGVector<DGstress>& S22, DGVector<DG>& vx_residual,
    DGVector<DG>& vy_residual);

/*!
 * LDG edge term <u_hat, tau.n> of the strain equation, multiplied with the local inverse
 * strain mass matrix and ADDED to E11/E12/E22.
 *   interior edges  : u_hat = (0.5+a) u1 + (0.5-a) u2
 *   Dirichlet edges : u_hat = 0, no contribution
 *   free edges      : u_hat = u, freeEdges[side] lists the cells with a free edge on that
 *                     side (0 bottom, 1 right, 2 top, 3 left), see freeEdgesOfMesh()
 */
template <int DG, int DGstress, int DGstrain>
void ldgVelocityFlux(const ParametricMesh& smesh,
    const ParametricDGMomentumMap<DG, DGstress, DGstrain>& pmap, const FloatType a,
    const std::array<std::vector<size_t>, 4>& freeEdges, const DGVector<DG>& vx,
    const DGVector<DG>& vy, DGVector<DGstrain>& E11, DGVector<DGstrain>& E12,
    DGVector<DGstrain>& E22);

/*!
 * The free edges: edges of ice cells on the outer mesh boundary which are not Dirichlet edges
 * (smesh.dirichlet). One list of cells per side (0 bottom, 1 right, 2 top, 3 left).
 */
std::array<std::vector<size_t>, 4> freeEdgesOfMesh(const ParametricMesh& smesh);

/*!
 * The LDG momentum equation with the mEVP solver.
 *
 * The Dirichlet edges are taken from smesh.dirichlet, which must be filled before, e.g. by
 * ParametricMesh::dirichletFromLandmask(false): u = 0 on edges between ice and land, the
 * outer mesh boundary is free (sigma.n = 0). CARTESIAN coordinates only, no periodic
 * boundaries.
 *
 * Usage per time step: prepareIteration(), then mEVPStep() nSteps times.
 */
template <int DG, int DGadvection> class LDGMomentum {
public:
    //! Precomputed matrices of the LDG discretization
    ParametricDGMomentumMap<DG, DG, DG> pmap;

    //! Ice velocity
    DGVector<DG> vx, vy;
    //! Ice velocity at the beginning of the time step
    DGVector<DG> vx0, vy0;
    //! Ocean and atmosphere (wind) velocity
    DGVector<DG> ox, oy, ax, ay;
    //! Strain and stress
    DGVector<DG> E11, E12, E22;
    DGVector<DG> S11, S12, S22;

    //! Gradient of the sea surface height in the GAUSSPOINTS(DG) Gauss points of each cell
    Eigen::Matrix<FloatType, Eigen::Dynamic, GAUSSPOINTS(DG), Eigen::RowMajor> sshGradX, sshGradY;

    //! Parameters of the PCG solver for the velocity
    FloatType pcgTolerance = 1.e-8; //!< relative residual
    int pcgMaxIterations = 1000;

    /*!
     * a and b are the parameters of the LDG fluxes: alternating flux u_hat = (0.5+a) u1 +
     * (0.5-a) u2 and penalty b/h. They are used in the jump penalty which is built once here.
     */
    LDGMomentum(const ParametricMesh& smesh, const FloatType a = 0.4, const FloatType b = 1.e9);

    //! Computes the gradient of the sea surface height (dG0) in the Gauss points, via CG1
    void setSeaSurfaceHeight(const DGVector<1>& ssh);

    /*!
     * Call once per time step, before the mEVP sub-iterations: stores the velocity of the
     * time step, the ice height (limited from below for the mass matrix) and concentration,
     * and sets up the velocity operator (1+beta) rho_ice M_H + dt J.
     */
    void prepareIteration(const VPParameters& params, const DGVector<DGadvection>& H,
        const DGVector<DGadvection>& A, const FloatType beta, const FloatType dt);

    /*!
     * One mEVP sub-iteration: (a) stress update, (b) velocity: PCG for
     * ((1+beta) rho_ice M_H + dt J) u = F, (c) strain. alpha, beta and dt must be the values
     * given to prepareIteration().
     */
    void mEVPStep(const VPParameters& params, const FloatType alpha, const FloatType beta,
        const FloatType dt);

    //! The strain of the current velocity: E = 0.5 (grad u + grad u^T), LDG
    void projectVelocityToStrain();

    //! Number of PCG iterations in the last velocity solve (max of x and y)
    int lastPCGIterations() const { return pcgIterations; }

private:
    const ParametricMesh& smesh;
    const FloatType a; //!< LDG flux parameter
    const std::array<std::vector<size_t>, 4> freeEdges; //!< see freeEdgesOfMesh()

    //! Ice height (limited from below) and concentration of the current time step
    DGVector<DGadvection> H, A;

    FloatType velocityMassFactor = 0.0; //!< (1+beta)*rho_ice
    FloatType velocityDt = 0.0; //!< dt multiplying the jump penalty J
    int pcgIterations = 0;
};

} /* namespace Nextsim */

#endif /* __LDGMOMENTUM_HPP */
