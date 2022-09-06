/*!
 * @file    ParametricTransport.hpp
 * @date    July 10, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICTRANSPORT_HPP
#define __PARAMETRICTRANSPORT_HPP

#include "cgVector.hpp"
#include "dgVector.hpp"

namespace Nextsim {

#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))

template <int DG>
void parametricTransportOperator(const ParametricMesh& smesh, const double dt, const DGVector<DG>& vx,
    const DGVector<DG>& vy, const EdgeVector<EDGEDOFS(DG)>& evx,
    const EdgeVector<EDGEDOFS(DG)>& evy, const DGVector<DG>& phi, DGVector<DG>& phiup);

/*!
 * Main class to manage the transport scheme on the parametric ParametricMesh
 *
 * template parameter DG, EDGEDOFS(DG) are number of local unknowns, that is
 * (1,2,3,6) on the cell and (1,2,3) on the edge
 */
template <int DG>
class ParametricTransport {
protected:
    //! spatial mesh.
    const ParametricMesh& smesh;

    //! reference to the current velocity
    DGVector<DG> velx, vely;

    //! Specifies the time stepping scheme [rk1, rk2, rk3]
    std::string timesteppingscheme;

    //! normal velocity in edges parallel to X- and Y-axis
    EdgeVector<EDGEDOFS(DG)> normalvel_X, normalvel_Y;

    //! temporary vectors for time stepping
    DGVector<DG> tmp1, tmp2, tmp3;

    //! Internal functions

    /*!
     * Performs one time step transporting phi with the Fwd-Euler Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk1(const double dt, DGVector<DG>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk2(const double dt, DGVector<DG>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk3(const double dt, DGVector<DG>& phi);

public:
    ParametricTransport(const ParametricMesh& mesh)
        : smesh(mesh)
        , timesteppingscheme("rk2")
    {
        if (!(smesh.nelements > 0)) {
            std::cerr << "ParametricTransport: The mesh must already be initialized!" << std::endl;
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
    }

    // Access members

    const DGVector<DG>& GetVx() const
    {
        return velx;
    }
    const DGVector<DG>& GetVy() const
    {
        return vely;
    }
    DGVector<DG>& GetVx()
    {
        return velx;
    }
    DGVector<DG>& GetVy()
    {
        return vely;
    }

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
    template <int CG>
    void prepareAdvection(const CGVector<CG>& cg_vx, const CGVector<CG>& cg_vy);

    /*!
     * Performs one time step transporting phi
     *
     * @params phi is the vector of values to be transported
     */
    void step(const double dt, DGVector<DG>& phi);
};

} /* namespace Nextsim */

#undef EDGEDOFS

#endif /* __PARAMETRICTRANSPORT_HPP */
