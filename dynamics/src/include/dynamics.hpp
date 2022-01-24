/*----------------------------   dynamics.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dynamics_H
#define __dynamics_H
/*----------------------------   dynamics.h     ---------------------------*/

#include "dgtransport.hpp"
#include "dgvector.hpp"
#include "mesh.hpp"

namespace Nextsim {
// pre-computed matrices for assembling dG-transport
// the Gauss rule for dG(n) is set to n+1 points
// this might not be enough?
#include "dgbasisfunctions_gausspoints.hpp"

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
public:
    Mesh mesh;

    /*!
   * Main variables for the ice model. 
   * ?? What are good DG spaces and combinations?
   */
    CellVector<2> vx, vy; //!< velocity fields
    CellVector<1> S11, S12, S22; //!< entries of (symmetric) stress tensor
    CellVector<1> E11, E12, E22; //!< entries of (symmetric) strain stress tensor

    CellVector<2> A, H; //!< ice height and ice concentration
    CellVector<0> D; //!< ice damage. ?? Really dG(0) ??

    CellVector<2> oceanX, oceanY; //!< ocean forcing. ?? Higher order??
    CellVector<0> atmX, atmY; //!< ocean forcing. ?? Higher order??

    //! temporary vectors for time stepping
    CellVector<2> tmpX, tmpY;

    /*! 
   * Subclasses for managing DG transport and the momentum problem
   */
    DGTransport<2> dgtransport;

    //public:
    /*! 
   * Constructor, initializes the subclasses and gives references to vectors
   * such as the velocity field
   */
    Dynamics()
        : dgtransport(vx, vy)
    {
    }

    //! Access
    Mesh& GetMesh() { return mesh; }
    CellVector<2>& GetTMPX() { return tmpX; }
    CellVector<2>& GetTMPY() { return tmpY; }

    CellVector<2>& GetVX() { return vx; }
    CellVector<2>& GetVY() { return vy; }
    CellVector<1>& GetS11() { return S11; }
    CellVector<1>& GetS12() { return S12; }
    CellVector<1>& GetS22() { return S22; }
    CellVector<1>& GetE11() { return E11; }
    CellVector<1>& GetE12() { return E12; }
    CellVector<1>& GetE22() { return E22; }
    CellVector<0>& GetD() { return D; }
    CellVector<2>& GetH() { return H; }
    CellVector<2>& GetA() { return A; }
    // TODO: First order need to implement stabilisations etc.
    //CellVector<1>& GetVX() { return vx; }
    //CellVector<1>& GetVY() { return vy; }
    //CellVector<0>& GetS11() { return S11; }
    //CellVector<0>& GetS12() { return S12; }
    //CellVector<0>& GetS21() { return S21; }
    //CellVector<0>& GetS22() { return S22; }
    //CellVector<0>& GetD() { return D; }
    //CellVector<1>& GetH() { return H; }
    //CellVector<1>& GetA() { return A; }

    CellVector<0>& GetAtmX() { return atmX; }
    CellVector<0>& GetAtmY() { return atmY; }
    CellVector<2>& GetOceanX() { return oceanX; }
    CellVector<2>& GetOceanY() { return oceanY; }

    /*! 
   * Sets important parameters, initializes these and that
   */
    void BasicInit();

    //////////////////////////////////////////////////

    /**!
   * controls the flow of the dynamical core
   */

    void momentumSymmetry();

    void advectionStep(const double dt);
    void momentumSubsteps(const double dt_momentum);
    void step(const double dt, const double dt_momentum);

    //! functions required to compute various terms of the momentum eq.

public:
    //! Computation of the stress tensor
    void addStressTensor(double scaleSigma);
    //! Computes the strain rate tensor 1/2(nabla v + nabla v^T) and stores in S11,S12,S22
    void computeStrainRateTensor();

private:
    void addStressTensorCell(double scaleSigma); //!< (S, \nabla Phi)
    void addStressTensorEdges(double scaleSigma); //!< < {{S}} , Phi >
    void addStressTensorBoundary(double scaleSigma); //!< same on boundary (from one side)
        //! consistency edge terms coming from integration by parts
    void stressTensorEdgesY(double scaleSigma, size_t c1, size_t c2)
    {
        const LocalEdgeVector<1> S11left( // left/right refers to the element
            S11(c1, 0) + 0.5 * S11(c1, 1),
            S11(c1, 2));
        const LocalEdgeVector<1> S11right(
            S11(c2, 0) - 0.5 * S11(c2, 1),
            S11(c2, 2));
        const LocalEdgeVector<1> S12left(
            S12(c1, 0) + 0.5 * S12(c1, 1),
            S12(c1, 2));
        const LocalEdgeVector<1> S12right(
            S12(c2, 0) - 0.5 * S12(c2, 1),
            S12(c2, 2));

        const LocalEdgeVector<2> avgS11 = 0.5 * (S11left + S11right) * BiGe13;
        const LocalEdgeVector<2> avgS12 = 0.5 * (S12left + S12right) * BiGe13;

        // {{S}} * [[Phi]]
        tmpX.block<1, 6>(c1, 0) -= scaleSigma / mesh.hx * avgS11 * BiG23_1;
        tmpX.block<1, 6>(c2, 0) += scaleSigma / mesh.hx * avgS11 * BiG23_3;
        tmpY.block<1, 6>(c1, 0) -= scaleSigma / mesh.hx * avgS12 * BiG23_1;
        tmpY.block<1, 6>(c2, 0) += scaleSigma / mesh.hx * avgS12 * BiG23_3;
    }
    void stressTensorEdgesX(double scaleSigma, size_t c1, size_t c2)
    {
        const LocalEdgeVector<1> S12bot( // bot/top refers to the element
            S12(c1, 0) + 0.5 * S12(c1, 2),
            S12(c1, 1));
        const LocalEdgeVector<1> S12top(
            S12(c2, 0) - 0.5 * S12(c2, 2),
            S12(c2, 1));
        const LocalEdgeVector<1> S22bot(
            S22(c1, 0) + 0.5 * S22(c1, 2),
            S22(c1, 1));
        const LocalEdgeVector<1> S22top(
            S22(c2, 0) - 0.5 * S22(c2, 2),
            S22(c2, 1));

        const LocalEdgeVector<2> avgS12 = 0.5 * (S12top + S12bot) * BiGe13;
        const LocalEdgeVector<2> avgS22 = 0.5 * (S22top + S22bot) * BiGe13;

        tmpX.block<1, 6>(c1, 0) -= scaleSigma / mesh.hy * avgS12 * BiG23_2;
        tmpX.block<1, 6>(c2, 0) += scaleSigma / mesh.hy * avgS12 * BiG23_0;
        tmpY.block<1, 6>(c1, 0) -= scaleSigma / mesh.hy * avgS22 * BiG23_2;
        tmpY.block<1, 6>(c2, 0) += scaleSigma / mesh.hy * avgS22 * BiG23_0;
    }

    //! consistency terms coming from integration by parts on boundaries
    void stressTensorBoundaryLeft(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<1> S11left( // left edge of the element
            S11(c, 0) - 0.5 * S11(c, 1),
            S11(c, 2));
        const LocalEdgeVector<1> S12left(
            S12(c, 0) - 0.5 * S12(c, 1),
            S12(c, 2));

        tmpX.block<1, 6>(c, 0) += scaleSigma / mesh.hx * S11left * BiGe13 * BiG23_3;
        tmpY.block<1, 6>(c, 0) += scaleSigma / mesh.hx * S12left * BiGe13 * BiG23_3;
    }
    void stressTensorBoundaryRight(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<1> S11right(
            S11(c, 0) + 0.5 * S11(c, 1),
            S11(c, 2));
        const LocalEdgeVector<1> S12right(
            S12(c, 0) + 0.5 * S12(c, 1),
            S12(c, 2));

        tmpX.block<1, 6>(c, 0) -= scaleSigma / mesh.hx * S11right * BiGe13 * BiG23_1;
        tmpY.block<1, 6>(c, 0) -= scaleSigma / mesh.hx * S12right * BiGe13 * BiG23_1;
    }
    void stressTensorBoundaryUpper(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<1> S12top(
            S12(c, 0) + 0.5 * S12(c, 2),
            S12(c, 1));
        const LocalEdgeVector<1> S22top(
            S22(c, 0) + 0.5 * S22(c, 2),
            S22(c, 1));

        tmpX.block<1, 6>(c, 0) -= scaleSigma / mesh.hy * S12top * BiGe13 * BiG23_2;
        tmpY.block<1, 6>(c, 0) -= scaleSigma / mesh.hy * S22top * BiGe13 * BiG23_2;
    }
    void stressTensorBoundaryLower(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<1> S12bot(
            S12(c, 0) - 0.5 * S12(c, 2),
            S12(c, 1));
        const LocalEdgeVector<1> S22bot(
            S22(c, 0) - 0.5 * S22(c, 2),
            S22(c, 1));

        tmpX.block<1, 6>(c, 0) += scaleSigma / mesh.hy * S12bot * BiGe13 * BiG23_0;
        tmpY.block<1, 6>(c, 0) += scaleSigma / mesh.hy * S22bot * BiGe13 * BiG23_0;
    }

    //! enforce continuity of the velocity by Nitsche
public:
    void velocityContinuity(double gamma); //!< gamma/h < [v], [Phi] > to enforce weak continuity of velocities
private:
    void velocityContinuityY(size_t c1, size_t c2, double gamma)
    {
        const LocalEdgeVector<2> leftX(vx(c1, 0) + 0.5 * vx(c1, 1) + 1. / 6. * vx(c1, 3),
            vx(c1, 2) + 0.5 * vx(c1, 5),
            vx(c1, 4));
        const LocalEdgeVector<2> leftY(vy(c1, 0) + 0.5 * vy(c1, 1) + 1. / 6. * vy(c1, 3),
            vy(c1, 2) + 0.5 * vy(c1, 5),
            vy(c1, 4));

        const LocalEdgeVector<2> rightX(vx(c2, 0) - 0.5 * vx(c2, 1) + 1. / 6. * vx(c2, 3),
            vx(c2, 2) - 0.5 * vx(c2, 5),
            vx(c2, 4));
        const LocalEdgeVector<2> rightY(vy(c2, 0) - 0.5 * vy(c2, 1) + 1. / 6. * vy(c2, 3),
            vy(c2, 2) - 0.5 * vy(c2, 5),
            vy(c2, 4));

        const LocalEdgeVector<2> jumpX = (rightX - leftX) * BiGe23;
        tmpX.block<1, 6>(c1, 0) += gamma / mesh.hx / mesh.hx * jumpX * BiG23_1;
        tmpX.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * jumpX * BiG23_3;

        const LocalEdgeVector<2> jumpY = (rightY - leftY) * BiGe23;

        tmpY.block<1, 6>(c1, 0) += gamma / mesh.hx / mesh.hx * jumpY * BiG23_1;
        tmpY.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * jumpY * BiG23_3;
    }
    void velocityContinuityX(size_t c1, size_t c2, double gamma)
    {
        const LocalEdgeVector<2> topX(
            vx(c1, 0) + 0.5 * vx(c1, 2) + 1. / 6. * vx(c1, 4),
            vx(c1, 1) + 0.5 * vx(c1, 5),
            vx(c1, 3));
        const LocalEdgeVector<2> bottomX(
            vx(c2, 0) - 0.5 * vx(c2, 2) + 1. / 6. * vx(c2, 4),
            vx(c2, 1) - 0.5 * vx(c2, 5), vx(c2, 3));

        const LocalEdgeVector<2> topY(
            vy(c1, 0) + 0.5 * vy(c1, 2) + 1. / 6. * vy(c1, 4),
            vy(c1, 1) + 0.5 * vy(c1, 5),
            vy(c1, 3));
        const LocalEdgeVector<2> bottomY(
            vy(c2, 0) - 0.5 * vy(c2, 2) + 1. / 6. * vy(c2, 4),
            vy(c2, 1) - 0.5 * vy(c2, 5), vy(c2, 3));

        const LocalEdgeVector<2> jumpX = (topX - bottomX) * BiGe23;
        tmpX.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * jumpX * BiG23_2;
        tmpX.block<1, 6>(c2, 0) += gamma / mesh.hy / mesh.hy * jumpX * BiG23_0;

        const LocalEdgeVector<2> jumpY = (topY - bottomY) * BiGe23;
        tmpY.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * jumpY * BiG23_2;
        tmpY.block<1, 6>(c2, 0) += gamma / mesh.hy / mesh.hy * jumpY * BiG23_0;
    }

    //! weak Dirichlet boundary with Nitsche method
public:
    void velocityDirichletBoundary(double gamma);

private:
    void velocityDirichletBoundaryLeft(const size_t c2, double gamma)
    { //x=0, y=t
        const LocalEdgeVector<2> rightX(vx(c2, 0) - 0.5 * vx(c2, 1) + 1. / 6. * vx(c2, 3),
            vx(c2, 2) - 0.5 * vx(c2, 5),
            vx(c2, 4));
        const LocalEdgeVector<2> rightY(vy(c2, 0) - 0.5 * vy(c2, 1) + 1. / 6. * vy(c2, 3),
            vy(c2, 2) - 0.5 * vy(c2, 5),
            vy(c2, 4));

        tmpX.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * rightX * BiGe23 * BiG23_3;
        tmpY.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * rightY * BiGe23 * BiG23_3;
    }

    void velocityDirichletBoundaryRight(const size_t c1, double gamma)
    { //x=1, y=t
        const LocalEdgeVector<2> leftX(vx(c1, 0) + 0.5 * vx(c1, 1) + 1. / 6. * vx(c1, 3),
            vx(c1, 2) + 0.5 * vx(c1, 5),
            vx(c1, 4));
        const LocalEdgeVector<2> leftY(vy(c1, 0) + 0.5 * vy(c1, 1) + 1. / 6. * vy(c1, 3),
            vy(c1, 2) + 0.5 * vy(c1, 5),
            vy(c1, 4));

        tmpX.block<1, 6>(c1, 0) -= gamma / mesh.hx / mesh.hx * leftX * BiGe23 * BiG23_1;
        tmpY.block<1, 6>(c1, 0) -= gamma / mesh.hx / mesh.hx * leftY * BiGe23 * BiG23_1;
    }

    void velocityDirichletBoundaryTop(const size_t c1, double gamma)
    {
        const LocalEdgeVector<2> topX(
            vx(c1, 0) + 0.5 * vx(c1, 2) + 1. / 6. * vx(c1, 4),
            vx(c1, 1) + 0.5 * vx(c1, 5),
            vx(c1, 3));

        const LocalEdgeVector<2> topY(
            vy(c1, 0) + 0.5 * vy(c1, 2) + 1. / 6. * vy(c1, 4),
            vy(c1, 1) + 0.5 * vy(c1, 5),
            vy(c1, 3));

        tmpX.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * topX * BiGe23 * BiG23_2;
        tmpY.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * topY * BiGe23 * BiG23_2;
    }

    void velocityDirichletBoundaryBottom(const size_t c2, double gamma)
    {
        const LocalEdgeVector<2> bottomX(
            vx(c2, 0) - 0.5 * vx(c2, 2) + 1. / 6. * vx(c2, 4),
            vx(c2, 1) - 0.5 * vx(c2, 5), vx(c2, 3));

        const LocalEdgeVector<2> bottomY(
            vy(c2, 0) - 0.5 * vy(c2, 2) + 1. / 6. * vy(c2, 4),
            vy(c2, 1) - 0.5 * vy(c2, 5), vy(c2, 3));

        tmpX.block<1, 6>(c2, 0) -= gamma / mesh.hy / mesh.hy * bottomX * BiGe23 * BiG23_0;
        tmpY.block<1, 6>(c2, 0) -= gamma / mesh.hy / mesh.hy * bottomY * BiGe23 * BiG23_0;
    }
};

}

/*----------------------------   dynamics.h     ---------------------------*/
/* end of #ifndef __dynamics_H */
#endif
/*----------------------------   dynamics.h     ---------------------------*/
