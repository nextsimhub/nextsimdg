/*----------------------------   dynamics.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dynamics_H
#define __dynamics_H
/*----------------------------   dynamics.h     ---------------------------*/

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

public:
    Dynamics()
    {
    }

    //! Performs one macro timestep1
    void macrostep()
    {
        advection_step();

        momentum_substeps();
    }
};

}

/*----------------------------   dynamics.h     ---------------------------*/
/* end of #ifndef __dynamics_H */
#endif
/*----------------------------   dynamics.h     ---------------------------*/
