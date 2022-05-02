/*!
 * @file    dgTransport.hpp
 * @date    Oct. 5, 2021
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGTRANSPORT_HPP
#define __DGTRANSPORT_HPP

#include "dgVector.hpp"

namespace Nextsim {

/*!
 * Main class to manage the DGdegree transport scheme
 *
 * template parameter DGdegree is degree (0,1,2) of dg scheme
 */
  template <int DGcell, int DGedge>
class DGTransport {
protected:
    //! reference to the current velocity
    const CellVector<DGcell>& velx;
    const CellVector<DGcell>& vely;

    //! Specifies the time stepping scheme [rk1, rk2, rk3]
    std::string timesteppingscheme;

    //! spatial mesh.
    Mesh mesh;

    //! velocity in edges
    EdgeVector<DGedge> velx_edgeY, vely_edgeX;

    //! temporary vectors for time stepping
    CellVector<DGcell> tmp1, tmp2, tmp3;

public:
    DGTransport(const CellVector<DGcell>& vx, const CellVector<DGcell>& vy)
        : velx(vx)
        , vely(vy)
        , timesteppingscheme("rk2")
    {
    }

    // High level functions

    void settimesteppingscheme(const std::string tss)
    {
        timesteppingscheme = tss;
        assert((tss == "rk1") || (tss == "rk2") || (tss == "rk3"));
    }

    /*!
     * Sets the spacial mesh
     *
     * Since the mesh has nearly no data, the entries are just copied.
     * The local vectors are initialized
     * - tmp1, tmp2, tmp3
     * - veledgex, veledgey
     */
    void setmesh(const Mesh& _mesh);

    /*!
     * Sets the velocity vector and prepares it for use in edges
     *
     * Stores pointer (?) to the velocity vector and creates the
     * corresponding structure on the edges, eg. velx and vely
     * and it also computes the zeros of the normal velocities
     * to be stored in velzerosx, velzerosy (normed to [0,1])
     *
     */
    void reinitvelocity();

    /*!
     * Performs one time step transporting phi with the Fwd-Euler Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk1(const double dt, CellVector<DGcell>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk2(const double dt, CellVector<DGcell>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk3(const double dt, CellVector<DGcell>& phi);

    /*!
     * Performs one time step transporting phi
     *
     * @params phi is the vector of values to be transported
     */
    void step(const double dt, CellVector<DGcell>& phi);
};

} /* namespace Nextsim */

#endif /* __DGTRANSPORT_HPP */
