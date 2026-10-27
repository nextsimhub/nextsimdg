/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGTRANSPORT_HPP
#define __DGTRANSPORT_HPP

#include "ParametricMap.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"

#include <algorithm>
#include <vector>

namespace Nextsim {

template <int DG>
void DGTransportOperator(const ParametricMesh& smesh, const FloatType dt, const DGVector<DG>& vx,
    const DGVector<DG>& vy, const EdgeVector<EDGEDOFS(DG)>& evx,
    const EdgeVector<EDGEDOFS(DG)>& evy, const DGVector<DG>& phi, DGVector<DG>& phiup);

/*!
 * Main class to manage the transport scheme on the parametric ParametricMesh
 *
 * template parameter DG, EDGEDOFS(DG) are number of local unknowns, that is
 * (1,2,3,6) on the cell and (1,2,3) on the edge
 */
template <int DG> class DGTransport {
protected:
    //! spatial mesh.
    const ParametricMesh& smesh;

    //! Precomputed stencil-like matrices for efficient numerical quadrature
    ParametricTransportMap<DG> parammap;

    //! reference to the current velocity
    DGVector<DG> velx, vely;

    //! Specifies the time stepping scheme [rk1, rk2, rk3]
    std::string timesteppingscheme;

    //! normal velocity in edges parallel to X- and Y-axis
    EdgeVector<EDGEDOFS(DG)> normalvel_X, normalvel_Y;

    /*!
     * true: the velocity is projected onto the advection space DG (prepareAdvection, CG
     * velocity): cell_term, edge_term_X/Y, boundary_lower/upper/left/right.
     * false: the full momentum velocity is used, no projection (prepareAdvectionDG, DG
     * velocity of the LDG momentum): cell_term_full, edge_term_X/Y_full, boundary_term_full
     */
    bool projectVelocity = true;

    //! Exterior value of phi on outer boundary edges with inflow (full velocity only)
    FloatType boundaryInflowValue = 0.0;

    //! Full velocity: Gauss points per edge, phi * (u.n) * psi has degree 2*deg(DG) + deg(u),
    //! at most 3 points
    static constexpr int NGP_EDGE = std::min(NGAUSSFORDEGREE1D(2 * DGDEGREE(DG) + 2), 3);

    //! Full velocity: Gauss points per direction in the cell, at most 3
    static constexpr int NGP_CELL = std::min(NGAUSSFORDEGREE1D(2 * DGDEGREE(DG) + 3), 3);

    //! (u_hat.n) * (edge length) in the NGP_EDGE Gauss points of each edge (one row per edge,
    //! numbering as normalvel_X/Y), normal pointing to +y (X-edges) or +x (Y-edges).
    //! Interior: u_hat = (0.5+a) u_c1 + (0.5-a) u_c2 (c1 = bottom/left cell), outer boundary:
    //! u_hat = trace of the cell velocity, 0 on edges touching land.
    Eigen::Matrix<FloatType, Eigen::Dynamic, NGP_EDGE, Eigen::RowMajor> edgeNormalVel_X,
        edgeNormalVel_Y;

    //! w_g * [J dT^{-T} grad psi]_x / _y in the NGP_CELL^2 Gauss points of each cell
    using CellTermMatrix = Eigen::Matrix<FloatType, DG, NGP_CELL * NGP_CELL>;
    std::vector<CellTermMatrix, Eigen::aligned_allocator<CellTermMatrix>> cellTermFullX,
        cellTermFullY;

    //! full velocity in the NGP_CELL^2 Gauss points of each cell (one row per cell)
    Eigen::Matrix<FloatType, Eigen::Dynamic, NGP_CELL * NGP_CELL, Eigen::RowMajor> cellVel_X,
        cellVel_Y;

    //! temporary vectors for time stepping
    DGVector<DG> tmp1, tmp2, tmp3;

    //! Internal functions

    /*!
     * Performs one time step transporting phi with the Fwd-Euler Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk1(const FloatType dt, DGVector<DG>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk2(const FloatType dt, DGVector<DG>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk3(const FloatType dt, DGVector<DG>& phi);

public:
    DGTransport(const ParametricMesh& mesh)
        : smesh(mesh)
        , parammap(mesh)
        , timesteppingscheme("rk2")
    {
        if (!(smesh.nelements > 0)) {
            std::cerr << "DGTransport: The mesh must already be initialized!" << std::endl;
            abort();
        }
        // Resize DG vectors for storing the velocity
        velx.resize_by_mesh(smesh);
        vely.resize_by_mesh(smesh);

        // resize tmp-vectors for time stepping
        tmp1.resize_by_mesh(smesh);
        tmp2.resize_by_mesh(smesh);
        tmp3.resize_by_mesh(smesh);

        // resize vectors to store the normal-velocity on the edges
        normalvel_Y.resize_by_mesh(smesh, EdgeType::Y);
        normalvel_X.resize_by_mesh(smesh, EdgeType::X);

        // initialize the mapping and set up required matrices
        parammap.InitializeAdvectionCellTerms();
        parammap.InitializeInverseDGMassMatrix();
    }

    // Set the DG velocity fields directly for testing
    void setVelocities(const DGVector<DG>& vx, const DGVector<DG>& vy)
    {
        velx = vx;
        vely = vy;
    }

private:
    const DGVector<DG>& GetVx() const { return velx; }
    const DGVector<DG>& GetVy() const { return vely; }
    DGVector<DG>& GetVx() { return velx; }
    DGVector<DG>& GetVy() { return vely; }

public:
    // High level functions
    void settimesteppingscheme(const std::string tss)
    {
        timesteppingscheme = tss;
        assert((tss == "rk1") || (tss == "rk2") || (tss == "rk3"));
    }

    /*!
     * Sets the normal-velocity vector on the edges
     * The normal velocity is scaled with the length of the edge,
     * this already serves as the integraiton weight
     */
    void reinitnormalvelocity();

    /*!
     * Prepares the advection step:
     * - interpolates CG velocity to DG
     * - initializes normal velocity on the edges
     */
    template <int CG> void prepareAdvection(const CGVector<CG>& cg_vx, const CGVector<CG>& cg_vy);

    /*!
     * Prepares the advection step with a DG velocity (LDG momentum), DGm components:
     * - cell term: the full velocity in the Gauss points
     * - interior edges: u_hat = (0.5+a) u_c1 + (0.5-a) u_c2 from the traces of the full
     *   velocity (a = 0: average)
     * - outer boundary edges: trace of the cell velocity, upwind with the exterior state
     *   boundaryInflowValue
     * - edges touching land: no flux
     * Periodic boundaries are not supported.
     */
    template <int DGm>
    void prepareAdvectionDG(
        const DGVector<DGm>& dg_vx, const DGVector<DGm>& dg_vy, const FloatType a = 0.0);

    /*!
     * Performs one time step transporting phi
     *
     * @params phi is the vector of values to be transported
     */
    void step(const FloatType dt, DGVector<DG>& phi);

private:
    /*!
     * Several internal functions
     */

    //! computes all integrals on the elements and the edges
    void DGTransportOperator(const ParametricMesh& smesh, const FloatType dt,
        const DGVector<DG>& vx, const DGVector<DG>& vy, const EdgeVector<EDGEDOFS(DG)>& normalvel_X,
        const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const DGVector<DG>& phi, DGVector<DG>& phiup);

    void edge_term_X(const ParametricMesh& smesh, const FloatType dt, DGVector<DG>& phiup,
        const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c1,
        const size_t c2, const size_t ie);
    void edge_term_Y(const ParametricMesh& smesh, const FloatType dt, DGVector<DG>& phiup,
        const DGVector<DG>& phi, const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c1,
        const size_t c2, const size_t ie);

    void cell_term(const ParametricMesh& smesh, FloatType dt, DGVector<DG>& phiup,
        const DGVector<DG>& phi, const DGVector<DG>& vx, const DGVector<DG>& vy, const size_t ic);

    //! upwind edge terms with the velocity edgeNormalVel_X/Y in the Gauss points
    void edge_term_X_full(DGVector<DG>& phiup, const DGVector<DG>& phi, const FloatType dt,
        const size_t c1, const size_t c2, const size_t ie);
    void edge_term_Y_full(DGVector<DG>& phiup, const DGVector<DG>& phi, const FloatType dt,
        const size_t c1, const size_t c2, const size_t ie);

    //! upwind term on the outer boundary edge on side Side (0 bottom, 1 right, 2 top, 3 left)
    //! of cell, the missing neighbour is replaced by boundaryInflowValue
    template <int Side>
    void boundary_term_full(
        DGVector<DG>& phiup, const DGVector<DG>& phi, const FloatType dt, const size_t cell);

    //! fills cellTermFullX/Y (once, in the first call of prepareAdvectionDG)
    void InitializeCellTermsFull();

    //! cell term with the full velocity cellVel_X/Y in the NGP_CELL^2 Gauss points
    void cell_term_full(
        DGVector<DG>& phiup, const DGVector<DG>& phi, const FloatType dt, const size_t eid);
};

} /* namespace Nextsim */

#endif /* __DGTRANSPORT_HPP */
