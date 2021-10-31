#include "dynamics.hpp"
#include "stopwatch.hpp"

namespace Nextsim
{
  extern Nextsim::Timer GlobalTimer;


      /*! 
   * Sets important parameters, initializes these and that
   */
  void Dynamics::BasicInit()
    {
        //! Degree of the time stepping in advection
        dgtransport.settimesteppingscheme("rk3");

        //! set meshes
        dgtransport.setmesh(mesh);
        dgtransport.settimemesh(timemesh);

        //! Init Vectors
        vx.resize_by_mesh(mesh);
        vy.resize_by_mesh(mesh);
        A.resize_by_mesh(mesh);
        H.resize_by_mesh(mesh);
        S11.resize_by_mesh(mesh);
        S12.resize_by_mesh(mesh);
        S22.resize_by_mesh(mesh);
        D.resize_by_mesh(mesh);

        oceanX.resize_by_mesh(mesh);
        oceanY.resize_by_mesh(mesh);
        atmX.resize_by_mesh(mesh);
        atmY.resize_by_mesh(mesh);
    }

    //////////////////////////////////////////////////

    /**!
   * controls the flow of the dynamical core
   */

  void Dynamics::advection_step()
    {
        GlobalTimer.start("dyn -- adv -- reinit");
        dgtransport.reinitvelocity();
        GlobalTimer.stop("dyn -- adv -- reinit");

        GlobalTimer.start("dyn -- adv -- step");
        dgtransport.step(A); // performs one time step with the 2nd Order Heun scheme
        dgtransport.step(H); // performs one time step with the 2nd Order Heun scheme
        GlobalTimer.stop("dyn -- adv -- step");
    }

  void Dynamics::momentum_substeps()
  {
    // vx.col(0) = 0.5 * oceanX.col(0) + 0.5 * atmX.col(0);
    // vy.col(0) = 0.5 * oceanY.col(0) + 0.5 * atmY.col(0);
  }
  
  //! Performs one macro timestep1
  void Dynamics::step()
  {
    GlobalTimer.start("dyn");
    GlobalTimer.start("dyn -- adv");
    advection_step();
    GlobalTimer.stop("dyn -- adv");
    
    GlobalTimer.start("dyn -- mom");
    momentum_substeps();
    GlobalTimer.stop("dyn -- mom");
    GlobalTimer.stop("dyn");
  }
  
  
}
