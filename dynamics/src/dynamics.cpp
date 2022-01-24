#include "dynamics.hpp"
#include "../../applications/benchmark_data.hpp"
#include "dgvisu.hpp"
#include "stopwatch.hpp"

namespace Nextsim {
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

    //! Init Vectors
    vx.resize_by_mesh(mesh);
    vy.resize_by_mesh(mesh);
    A.resize_by_mesh(mesh);
    H.resize_by_mesh(mesh);
    S11.resize_by_mesh(mesh);
    S12.resize_by_mesh(mesh);
    S22.resize_by_mesh(mesh);

    E11.resize_by_mesh(mesh);
    E12.resize_by_mesh(mesh);
    E22.resize_by_mesh(mesh);
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

void Dynamics::advectionStep(const double dt)
{
    GlobalTimer.start("dyn -- adv -- reinit");
    dgtransport.reinitvelocity();
    GlobalTimer.stop("dyn -- adv -- reinit");

    GlobalTimer.start("dyn -- adv -- step");
    dgtransport.step(dt, A); // performs one time step with the 2nd Order Heun scheme
    dgtransport.step(dt, H); // performs one time step with the 2nd Order Heun scheme

    // Limit H and A
#pragma omp parallel for
    for (size_t i = 0; i < mesh.n; ++i) {
        A(i, 0) = std::max(std::min(A(i, 0), 1.0), 0.0);
        H(i, 0) = std::max(H(i, 0), 0.0);
    }

    GlobalTimer.stop("dyn -- adv -- step");
}

void Dynamics::velocityContinuity(double gamma)
{
    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            velocityContinuityY(ic, ic + 1, gamma);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            velocityContinuityX(ic, ic + mesh.nx, gamma);
    }
}

void Dynamics::addStressTensorEdges(double scaleSigma)
{

    // Y - edges, only inner ones
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        size_t ic = iy * mesh.nx; // first index of left cell in row

        for (size_t i = 0; i < mesh.nx - 1; ++i, ++ic)
            stressTensorEdgesY(scaleSigma, ic, ic + 1);
    }

    // X - edges, only inner ones
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {
        size_t ic = ix; // first index of left cell in column
        for (size_t i = 0; i < mesh.ny - 1; ++i, ic += mesh.nx)
            stressTensorEdgesX(scaleSigma, ic, ic + mesh.nx);
    }
}

void Dynamics::addStressTensorBoundary(double scaleSigma)
{
    //! consistency term on lower and upper boundary
#pragma omp parallel for
    for (size_t ix = 0; ix < mesh.nx; ++ix) {

        const size_t clower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;

        stressTensorBoundaryUpper(scaleSigma, cupper);
        stressTensorBoundaryLower(scaleSigma, clower);
    }

    //! consistency term on left and right boundary
#pragma omp parallel for
    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t cright = (iy + 1) * mesh.nx - 1;

        stressTensorBoundaryLeft(scaleSigma, cleft);
        stressTensorBoundaryRight(scaleSigma, cright);
    }
}

void Dynamics::momentumSymmetry() { }

void Dynamics::velocityDirichletBoundary(double gamma)
{

    for (size_t ix = 0; ix < mesh.nx; ++ix) {

        const size_t clower = ix;
        const size_t cupper = mesh.n - mesh.nx + ix;

        velocityDirichletBoundaryTop(cupper, gamma);
        velocityDirichletBoundaryBottom(clower, gamma);
    }

    for (size_t iy = 0; iy < mesh.ny; ++iy) {
        const size_t cleft = iy * mesh.nx;
        const size_t cright = (iy + 1) * mesh.nx - 1;

        velocityDirichletBoundaryLeft(cleft, gamma);
        velocityDirichletBoundaryRight(cright, gamma);
    }
}

void Dynamics::computeStrainRateTensor()
{
    E11.col(0) = 1. / mesh.hx * vx.col(1);
    E11.col(1) = 1. / mesh.hx * 2. * vx.col(3);
    E11.col(2) = 1. / mesh.hx * vx.col(5);
    E12.col(0) = 1. * 0.5 * (vy.col(1) / mesh.hx + vx.col(2) / mesh.hy);
    E12.col(1) = 1. * 0.5 * (vx.col(5) / mesh.hy + 2.0 * vy.col(3) / mesh.hx);
    E12.col(2) = 1. * 0.5 * (2.0 * vx.col(4) / mesh.hy + vy.col(5) / mesh.hx);
    E22.col(0) = 1. / mesh.hy * vy.col(2);
    E22.col(1) = 1. / mesh.hy * vy.col(5);
    E22.col(2) = 1. / mesh.hy * 2. * vy.col(4);
}

//! Computes the cell terms of sigma = 1/2(nabla v + nabla v^T)
void Dynamics::addStressTensorCell(double scaleSigma)
{
    // S11 d_x phi_x + S12 d_y phi_x
    tmpX.col(1) += 12. * scaleSigma / mesh.hx * (S11.col(0));
    tmpX.col(3) += 180. * scaleSigma / mesh.hx * (S11.col(1) / 6.);
    tmpX.col(5) += 144. * scaleSigma / mesh.hx * (S11.col(2) / 12.);

    tmpX.col(2) += 12. * scaleSigma / mesh.hy * (S12.col(0));
    tmpX.col(4) += 180. * scaleSigma / mesh.hy * (S12.col(2) / 6.);
    tmpX.col(5) += 144. * scaleSigma / mesh.hy * (S12.col(1) / 12.);

    // S12 d_x phi_y + S22 d_y phi_y
    tmpY.col(1) += 12. * scaleSigma / mesh.hx * (S12.col(0));
    tmpY.col(3) += 180. * scaleSigma / mesh.hx * (S12.col(1) / 6.);
    tmpY.col(5) += 144. * scaleSigma / mesh.hx * (S12.col(2) / 12.);

    tmpY.col(2) += 12. * scaleSigma / mesh.hy * (S22.col(0));
    tmpY.col(4) += 180. * scaleSigma / mesh.hy * (S22.col(2) / 6.);
    tmpY.col(5) += 144. * scaleSigma / mesh.hy * (S22.col(1) / 12.);
}

void Dynamics::addStressTensor(double scaleSigma)
{
    addStressTensorCell(scaleSigma); //!< cell term
    addStressTensorEdges(scaleSigma);
    addStressTensorBoundary(scaleSigma);
}

void Dynamics::momentumSubsteps(const double dt_momentum)
{
    //MOMENTUM EQUATION
    tmpX.zero();
    tmpY.zero();

    // d_t U = ...

    // ocean
    // L/(rho H) * ReferenceScale::C_ocean * ReferenceScale::rho_ocean * |velwater - vel| (velwater-vel)
    tmpX.col(0) += 1 / ReferenceScale::rho_ice * ReferenceScale::C_ocean * ReferenceScale::rho_ocean * ((oceanX.col(0) - vx.col(0)).array().abs() / H.col(0).array() * oceanX.col(0).array()).matrix();

    tmpX -= 1 / ReferenceScale::rho_ice * ReferenceScale::C_ocean * ReferenceScale::rho_ocean * (vx.array().colwise() * ((oceanX.col(0) - vx.col(0)).array().abs() / H.col(0).array())).matrix();

    tmpY.col(0) += 1 / ReferenceScale::rho_ice * ReferenceScale::C_ocean * ReferenceScale::rho_ocean * ((oceanY.col(0) - vy.col(0)).array().abs() / H.col(0).array() * oceanY.col(0).array()).matrix();

    tmpY -= 1 / ReferenceScale::rho_ice * ReferenceScale::C_ocean * ReferenceScale::rho_ocean * (vy.array().colwise() * ((oceanY.col(0) - vy.col(0)).array().abs() / H.col(0).array())).matrix();

    // atm.
    // L/(rho H) * ReferenceScale::C_atm * ReferenceScale::rho_atm * |velatm| velatm
    tmpX.col(0) += 1 / ReferenceScale::rho_ice * ReferenceScale::C_atm * ReferenceScale::rho_atm * (atmX.col(0).array().abs() / H.col(0).array() * atmX.col(0).array()).matrix();
    tmpY.col(0) += 1 / ReferenceScale::rho_ice * ReferenceScale::C_atm * ReferenceScale::rho_atm * (atmY.col(0).array().abs() / H.col(0).array() * atmY.col(0).array()).matrix();

    /**! 
     *
     * Compute the stress tensor ...
     * ...
     *
     */
    const double scaleSigma = -1.; //!< -1 since -div(S)  term goes to rhs

    // This is already defined in benchmark.cpp
    const double gamma = 5.0; //!< parameter in front of internal penalty terms
    const double gammaboundary = 5.0; //!< parameter in front of boundary penalty terms

    addStressTensor(scaleSigma); //!<  +div(S, nabla phi) - < Sn, phi>
    velocityContinuity(gamma); //!< penalize velocity jump in inner edges
    velocityDirichletBoundary(gammaboundary); //!< no-slip for velocity

    vx += dt_momentum * tmpX;
    vy += dt_momentum * tmpY;

    //DAMAGE EQUATION

    // Sigma = D = sym(grad v)
    computeStrainRateTensor(); //!< S = 1/2 (nabla v + nabla v^T)

    GlobalTimer.start("dyn -- mom -- damage");
#pragma omp parallel for
    for (size_t i = 0; i < mesh.n; ++i) {
        // Compute Pmax Eqn.(8) the way like in nextsim finiteelement.cpp
        double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
        //std::cout << Pmax << " " << sigma_n << std::endl;
        double const expC = std::exp(ReferenceScale::compaction_param * (1. - A(i, 0)));
        double const time_viscous = ReferenceScale::undamaged_time_relaxation_sigma * std::pow((1. - D(i, 0)) * expC, ReferenceScale::exponent_relaxation_sigma - 1.);

        // Plastic failure tildeP
        double tildeP;
        if (sigma_n < 0.) {
            //below line copied from nextsim finiteelement.cpp
            //double const Pmax = std::pow(M_thick[cpt], exponent_compression_factor)*compression_factor*expC;
            double const Pmax = ReferenceScale::Pstar * pow(H(i, 0), 1.5) * exp(-20.0 * (1.0 - A(i, 0)));
            // tildeP must be capped at 1 to get an elastic response
            tildeP = std::min(1., -Pmax / sigma_n);
        } else {
            tildeP = 0.;
        }

        // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
        // min and - from nextsim
        double const multiplicator = std::min(1. - 1e-12,
            time_viscous / (time_viscous + dt_momentum * (1. - tildeP)));

        double const elasticity = ReferenceScale::young * (1. - D(i, 0)) * expC;

        double const Dunit_factor = 1. / (1. - SQR(ReferenceScale::nu0));
        /* Stiffness matrix
        * 1  nu 0
        * nu 1  0
        * 0  0  (1-nu)
        */
        //M_Dunit[0] = Dunit_factor * 1.;
        //M_Dunit[1] = Dunit_factor * nu0;
        //M_Dunit[3] = Dunit_factor * nu0;
        //M_Dunit[4] = Dunit_factor * 1.;
        //M_Dunit[8] = Dunit_factor * (1. - nu0) ;

        //compute SigmaE
        double SigmaE11 = Dunit_factor * (E11(i, 0) + ReferenceScale::nu0 * E22(i, 0));
        double SigmaE22 = Dunit_factor * (ReferenceScale::nu0 * E11(i, 0) + E22(i, 0));
        double SigmaE12 = 1. / (1 - ReferenceScale::nu0) * E12(i, 0);

        //Elasit prediction Eqn. (32)
        S11(i, 0) += dt_momentum * elasticity * SigmaE11;
        S11(i, 0) *= multiplicator;
        S12(i, 0) += dt_momentum * elasticity * SigmaE12;
        S12(i, 0) *= multiplicator;
        S22(i, 0) += dt_momentum * elasticity * SigmaE22;
        S22(i, 0) *= multiplicator;

        //continiue if stress in inside the failure envelope
        /* Compute the shear and normal stresses, which are two invariants of the internal stress tensor */
        double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
        //update sigma_n
        sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
        //cohesion Eqn. (21)
        //Reference length scale is fixed 0.1 since its cohesion parameter at the lab scale (10 cm)
        double const C_fix = ReferenceScale::C_lab * std::sqrt(0.1 / mesh.hx);
        ; // C_lab;...  : cohesion (Pa)

        // d critical Eqn. (29)
        double dcrit;
        if (sigma_n < -ReferenceScale::compr_strength)
            dcrit = -ReferenceScale::compr_strength / sigma_n;
        else
            // M_Cohesion[cpt] depends on local random contribution
            // M_Cohesion[i] = C_fix+C_alea*(M_random_number[i]);
            // finiteelement.cpp#L3834
            dcrit = C_fix / (sigma_s + ReferenceScale::tan_phi * sigma_n);

        /* Calculate the adjusted level of damage */
        if ((0. < dcrit) && (dcrit < 1.)) // sigma_s - tan_phi*sigma_n < 0 is always inside, but gives dcrit < 0
        {
            /* Calculate the characteristic time for damage and damage increment */
            // M_delta_x[cpt] = mesh.hx ???
            double const td = mesh.hx * std::sqrt(2. * (1. + ReferenceScale::nu0) * ReferenceScale::rho_ice)
                / std::sqrt(elasticity);

            // Eqn. (34)
            D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

            // Recalculate the new state of stress by relaxing elstically Eqn. (36)
            S11(i, 0) -= S11(i, 0) * (1. - dcrit) * dt_momentum / td;
            S12(i, 0) -= S12(i, 0) * (1. - dcrit) * dt_momentum / td;
            S22(i, 0) -= S22(i, 0) * (1. - dcrit) * dt_momentum / td;
        }
    }
    GlobalTimer.stop("dyn -- mom -- damage");
}

//! Performs one macro timestep1
void Dynamics::step(const double dt, const double dt_momentum)
{
    GlobalTimer.start("dyn");
    GlobalTimer.start("dyn -- adv");
    advectionStep(dt);
    GlobalTimer.stop("dyn -- adv");

    GlobalTimer.start("dyn -- mom");
    momentumSubsteps(dt_momentum);
    GlobalTimer.stop("dyn -- mom");
    GlobalTimer.stop("dyn");
}

}
