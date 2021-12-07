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

        tmpX.resize_by_mesh(mesh);
        tmpY.resize_by_mesh(mesh);

    }

    //////////////////////////////////////////////////

    /**!
   * controls the flow of the dynamical core
   */

  void Dynamics::advectionStep()
    {
        GlobalTimer.start("dyn -- adv -- reinit");
        dgtransport.reinitvelocity();
        GlobalTimer.stop("dyn -- adv -- reinit");

        GlobalTimer.start("dyn -- adv -- step");
        dgtransport.step(A); // performs one time step with the 2nd Order Heun scheme
        dgtransport.step(H); // performs one time step with the 2nd Order Heun scheme
        GlobalTimer.stop("dyn -- adv -- step");
    }

    void Dynamics::momentumJumps()
    {
          // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row
        
        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            stabilizationY(ic, ic + 1);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            stabilizationX(ic, ic + mesh.nx);
    }

    }

    void Dynamics::momentumDirichletBoundary()
    {
      
      for (size_t ix = 0; ix < mesh.nx; ++ix) {

          const size_t clower = ix;
          const size_t cupper = mesh.n - mesh.nx + ix;

          boundaryDirichletTop<2>(cupper);
          boundaryDirichletBottom<2>(clower);

      }

      for (size_t iy = 0; iy < mesh.ny; ++iy) {
          const size_t cleft = iy * mesh.nx;
          const size_t cright = (iy + 1) * mesh.nx - 1 ;

          boundaryDirichletLeft<2>(cleft);
          boundaryDirichletRight<2>(cright);

      }

    }//momentumBoundaryStabilization

  void Dynamics::momentumSubsteps()
  {
    tmpX.zero();
    tmpY.zero();

    double L = 512000.0;
    double T = 1000.0;
    double Cwater = 5.5e-3;
    double Catm = 1.2e-3;
    double rhoWater = 1026.0;
    double rhoIce = 900.0;
    double rhoAtm = 1.3;
    double Sigma  = 1;

    // d_t U = ...

// ocean

    // L/(rho H) * Cwater * rhoWater * |velwater - vel| (velwater-vel)
    tmpX.col(0) += L / rhoIce * Cwater * rhoWater *
     ((oceanX.col(0)-vx.col(0)).array().abs()/ 
       H.col(0).array() * 
       oceanX.col(0).array()).matrix();

   tmpX -= L / rhoIce * Cwater * rhoWater *
     (vx.array().colwise() * ((oceanX.col(0)-vx.col(0)).array().abs()/ H.col(0).array())).matrix();

    tmpY.col(0) += L / rhoIce * Cwater * rhoWater *
     ((oceanY.col(0)-vy.col(0)).array().abs()/ 
       H.col(0).array() * 
       oceanY.col(0).array()).matrix();

   tmpY -= L / rhoIce * Cwater * rhoWater *
     (vy.array().colwise() * ((oceanY.col(0)-vy.col(0)).array().abs()/ H.col(0).array())).matrix();


   // atm.
    // L/(rho H) * Catm * rhoAtm * |velatm| velatm
    tmpX.col(0) += L / rhoIce * Catm * rhoAtm *
     (atmX.col(0).array().abs()/ H.col(0).array()
      * atmX.col(0).array()).matrix();

    tmpY.col(0) += L / rhoIce * Catm * rhoAtm *
     (atmY.col(0).array().abs()/ H.col(0).array()
      * atmY.col(0).array()).matrix();

  const double scaleSigma = Sigma * T * T / L / L / rhoIce;
  //( sigma,grad phi ) DG0 component
  tmpX.col(0) += scaleSigma*S11.col(0)  ;
  tmpY.col(0) += scaleSigma*S22.col(0)  ;

  // jump stabilization
  momentumJumps();

  // boundary zero Dirichlet
  momentumDirichletBoundary();

   vx += timemesh.dt_momentum * tmpX;
   vy += timemesh.dt_momentum * tmpY;

    // vx.col(0) = 0.5 * oceanX.col(0) + 0.5 * atmX.col(0);
    // vy.col(0) = 0.5 * oceanY.col(0) + 0.5 * atmY.col(0);
  }
  
  //! Performs one macro timestep1
  void Dynamics::step()
  {
    GlobalTimer.start("dyn");
    GlobalTimer.start("dyn -- adv");
    advectionStep();
    GlobalTimer.stop("dyn -- adv");
    
    GlobalTimer.start("dyn -- mom");
    momentumSubsteps();
    GlobalTimer.stop("dyn -- mom");
    GlobalTimer.stop("dyn");
  }
  
  
}
