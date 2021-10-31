/*----------------------------   dynamics.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dynamics_H
#define __dynamics_H
/*----------------------------   dynamics.h     ---------------------------*/

#include "dgtransport.hpp"
#include "dgvector.hpp"
#include "mesh.hpp"
#include "timemesh.hpp"

namespace Nextsim {

/*!
   * This class controls the timestepping of the dynamical core
   * - Advection
   * - Subcycling of momentum and damage
   *
   * This class is the main class that contains all the required data
   */
class Dynamics {

    /*!
     * Here we define mesh and time mesh. References
     * to these classes will be used in the submodules
     */
    Mesh mesh;
    TimeMesh timemesh;

    /*!
   * Main variables for the ice model. 
   * ?? What are good DG spaces and combinations?
   */
    CellVector<2> vx, vy; //!< velocity fields
    CellVector<1> S11, S12, S22; //!< entries of (symmetric) stress tensor
    CellVector<2> A, H; //!< ice height and ice concentration
    CellVector<0> D; //!< ice damage. ?? Really dG(0) ??

    CellVector<0> oceanX, oceanY; //!< ocean forcing. ?? Higher order??
    CellVector<0> atmX, atmY; //!< ocean forcing. ?? Higher order??

    /*! 
   * Subclasses for managing DG transport and the momentum problem
   */
    DGTransport<2> dgtransport;

public:
    /*! 
   * Constructor, initializes the subclasses and gives references to vectors
   * such as the velocity field
   */
    Dynamics()
        : dgtransport(vx, vy)
    {
    }

    //! Access
    Mesh& GetMesh()
    {
        return mesh;
    }
    TimeMesh& GetTimeMesh()
    {
        return timemesh;
    }
    CellVector<2>& GetVX()
    {
        return vx;
    }
    CellVector<2>& GetVY()
    {
        return vy;
    }
    CellVector<1>& GetS11()
    {
        return S11;
    }
    CellVector<1>& GetS12()
    {
        return S12;
    }
    CellVector<1>& GetS22()
    {
        return S22;
    }
    CellVector<0>& GetD()
    {
        return D;
    }
    CellVector<2>& GetH()
    {
        return H;
    }
    CellVector<2>& GetA()
    {
        return A;
    }

    CellVector<0>& GetAtmX()
    {
        return atmX;
    }
    CellVector<0>& GetAtmY()
    {
        return atmY;
    }
    CellVector<0>& GetOceanX()
    {
        return oceanX;
    }
    CellVector<0>& GetOceanY()
    {
        return oceanY;
    }

    /*! 
   * Sets important parameters, initializes these and that
   */
    void BasicInit();

    //////////////////////////////////////////////////

    /**!
   * controls the flow of the dynamical core
   */

    void advection_step();
    void momentum_substeps();
    void step();
};

}

/*----------------------------   dynamics.h     ---------------------------*/
/* end of #ifndef __dynamics_H */
#endif
/*----------------------------   dynamics.h     ---------------------------*/
