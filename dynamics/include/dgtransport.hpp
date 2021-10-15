/*----------------------------   dgtransport.h     ---------------------------*/
#ifndef __dgtransport_H
#define __dgtransport_H
/*----------------------------   dgtransport.h     ---------------------------*/

/*!
 * @file    dgtransport.h
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @date    Oct. 5, 2021
 */

#include "dgvector.hpp"
#include "timemesh.hpp"

namespace Nextsim {

/*!
 * Main class to manage the DGdegree transport scheme
 *
 * template parameter DGdegree is degree (0,1,2) of dg scheme
 */
template <int DGdegree>
class DGTransport {
protected:
    //! pointer to the current velocity
    const CellVector<DGdegree>* velxpointer;
    const CellVector<DGdegree>* velypointer;

    //! spatial mesh.
    Mesh mesh;
    //! time mesh.
    TimeMesh timemesh;

    //! velocity in edges
    EdgeVector<DGdegree> xvel_edgeY, yvel_edgeX;

    //! temporary vectors for time stepping
    CellVector<DGdegree> tmp1, tmp2, tmp3;

public:
    DGTransport()
        : velxpointer(NULL)
        , velypointer(NULL)
    {
    }

    // High level functions

    /*!
   * Sets the spacial mesh
   *
   * Since the mesh has nearly no data, the entries are just copied.
   * The local vectors are initialized
   * - tmp1, tmp2, tmp3
   * - veledgex, veledgey
   */
    void setmesh(const Mesh& _mesh);

    //! Sets the timemesh. Just a copy
    void settimemesh(const TimeMesh& _timemesh);

    /*!
   * Sets the velocity vector and prepares it for use in edges
   *
   * Stores pointer (?) to the velocity vector and creates the
   * corresponding structure on the edges, eg. velx and vely
   * and it also computes the zeros of the normal velocities
   * to be stored in velzerosx, velzerosy (normed to [0,1])
   *
   */
    void setvelocity(const CellVector<DGdegree>& velx,
        const CellVector<DGdegree>& vely);

    /*!
   * Performs one time step transporting phi with the Fwd-Euler Scheme
   *
   * @params phi is the vector of values to be transported
   */
    void step_fwdeuler(CellVector<DGdegree>& phi);

    /*!
   * Performs one time step transporting phi with the 2nd Order Heun Scheme
   *
   * @params phi is the vector of values to be transported
   */
    void step_heun(CellVector<DGdegree>& phi);
};

} // namespace Nextsim

/*----------------------------   dgtransport.h     ---------------------------*/
/* end of #ifndef __dgtransport_H */
#endif
/*----------------------------   dgtransport.h     ---------------------------*/
