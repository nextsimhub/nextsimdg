/*----------------------------   dynamics.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __newdynamics_H
#define __newdynamics_H
/*----------------------------   dynamics.h     ---------------------------*/

#include "dgvector.hpp"

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
class NewDynamics {

public:
    void momentumSymmetry();

    void advectionStep();
    void momentumSubsteps();
    void step();

    //! functions required to compute various terms of the momentum eq.

public:
    //! Computation of the stress tensor
    void addStressTensor(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22);

    //! Computes the strain rate tensor 1/2(nabla v + nabla v^T) and stores in E11,E12,E22
    void computeStrainRateTensor(const Mesh& mesh, CellVector<2>& vx,
        CellVector<2>& vy, CellVector<1>& E11, CellVector<1>& E12, CellVector<1>& E22);

private:
    void addStressTensorCell(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22); //!< (S, \nabla Phi)
    void addStressTensorEdges(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22); //!< < {{S}} , Phi >
    void addStressTensorBoundary(double scaleSigma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22); //!< same on boundary (from one side)
    //! consistency edge terms coming from integration by parts

    void stressTensorEdgesY(double scaleSigma, size_t c1, size_t c2, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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
    void stressTensorEdgesX(double scaleSigma, size_t c1, size_t c2, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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
    void stressTensorBoundaryLeft(double scaleSigma, size_t c, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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
    void stressTensorBoundaryRight(double scaleSigma, size_t c, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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
    void stressTensorBoundaryUpper(double scaleSigma, size_t c, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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
    void stressTensorBoundaryLower(double scaleSigma, size_t c, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<1>& S11, CellVector<1>& S12, CellVector<1>& S22)
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
    void velocityContinuity(double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy); //!< gamma/h < [v], [Phi] > to enforce weak continuity of velocities
private:
    void velocityContinuityY(size_t c1, size_t c2, double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
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
    void velocityContinuityX(size_t c1, size_t c2, double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
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
    void velocityDirichletBoundary(double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy);

private:
    void velocityDirichletBoundaryLeft(const size_t c2, double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
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

    void velocityDirichletBoundaryRight(const size_t c1, double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
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

    void velocityDirichletBoundaryTop(const size_t c1, double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
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

    void velocityDirichletBoundaryBottom(const size_t c2, double gamma, const Mesh& mesh, CellVector<2>& tmpX,
        CellVector<2>& tmpY, CellVector<2>& vx, CellVector<2>& vy)
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
/* end of #ifndef __newdynamics_H */
#endif
/*----------------------------   dynamics.h     ---------------------------*/
