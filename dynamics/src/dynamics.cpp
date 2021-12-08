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


  void Dynamics::momentumConsistency() {

          // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row
        
        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            consistencyY(ic, ic + 1);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            consistencyX(ic, ic + mesh.nx);
  }

  }

  void Dynamics::momentumSymmetry() {}


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
    double Sigma  = 1e5;

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


  // sigma = D = sym(grad v)
  S11.col(0) = vx.col(1);
  S11.col(1) = 0.5*vx.col(3);
  S11.col(2) = vx.col(5);
  S12.col(0) = 0.5*(vy.col(1) + vx.col(2)) ;
  S12.col(1) = 0.5*( 0.5*vy.col(3) + vx.col(5) );
  S12.col(2) = 0.5*( vy.col(5) + 0.5*vx.col(4));
  S22.col(0) = vy.col(2);
  S22.col(1) = vy.col(5);
  S22.col(2) = 0.5*vy.col(4);

  const double scaleSigma = Sigma * T * T / L / L / rhoIce;
  //( sigma,grad phi ) 
  // S11 d_x phi_x + S12 d_y phi_x
  tmpX.col(1) += scaleSigma*(S11.col(0)    );
  tmpX.col(3) += scaleSigma*(S11.col(1)/6. );
  tmpX.col(5) += scaleSigma*(S11.col(2)/12.);
  tmpX.col(2) += scaleSigma*(S12.col(0)    );
  tmpX.col(4) += scaleSigma*(S12.col(2)/6. );
  tmpX.col(5) += scaleSigma*(S12.col(1)/12.);
  // S12 d_x phi_y + S22 d_y phi_y
  tmpY.col(1) += scaleSigma*(S12.col(0)    );
  tmpY.col(3) += scaleSigma*(S12.col(1)/6. );
  tmpY.col(5) += scaleSigma*(S12.col(2)/12.);
  tmpY.col(2) += scaleSigma*(S22.col(0)    );
  tmpY.col(4) += scaleSigma*(S22.col(2)/6. );
  tmpY.col(5) += scaleSigma*(S22.col(1)/12.);
 
  // Stress consistency and symmetry terms
  //avg(sigma) n jump(phi)
  momentumConsistency();
  //jump(v) n avg(grad phi)
  momentumSymmetry();



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
