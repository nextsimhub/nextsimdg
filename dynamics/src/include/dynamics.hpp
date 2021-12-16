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
    TimeMesh timemesh;

    /*!
   * Main variables for the ice model. 
   * ?? What are good DG spaces and combinations?
   */
    CellVector<2> vx, vy; //!< velocity fields
    CellVector<1> S11, S12, S22; //!< entries of (symmetric) stress tensor
    CellVector<1> S21; // for testing
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
    TimeMesh& GetTimeMesh() { return timemesh; }
    CellVector<2>& GetVX() { return vx; }
    CellVector<2>& GetVY() { return vy; }
    CellVector<1>& GetS11() { return S11; }
    CellVector<1>& GetS12() { return S12; }
    CellVector<1>& GetS21() { return S21; }
    CellVector<1>& GetS22() { return S22; }
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

    void stabilizationY(size_t c1, size_t c2)
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
        tmpX.block<1, 6>(c1, 0) += 1. / mesh.hx / mesh.hx * jumpX * BiG23_1;
        tmpX.block<1, 6>(c2, 0) -= 1. / mesh.hx / mesh.hx * jumpX * BiG23_3;

        const LocalEdgeVector<2> jumpY = (rightY - leftY) * BiGe23;

        tmpY.block<1, 6>(c1, 0) += 1. / mesh.hx / mesh.hx * jumpY * BiG23_1;
        tmpY.block<1, 6>(c2, 0) -= 1. / mesh.hx / mesh.hx * jumpY * BiG23_3;
    }

    void stabilizationX(size_t c1, size_t c2)
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
        tmpX.block<1, 6>(c1, 0) -= 1. / mesh.hy / mesh.hy * jumpX * BiG23_2;
        tmpX.block<1, 6>(c2, 0) += 1. / mesh.hy / mesh.hy * jumpX * BiG23_0;

        const LocalEdgeVector<2> jumpY = (topY - bottomY) * BiGe23;
        tmpY.block<1, 6>(c1, 0) -= 1. / mesh.hy / mesh.hy * jumpY * BiG23_2;
        tmpY.block<1, 6>(c2, 0) += 1. / mesh.hy / mesh.hy * jumpY * BiG23_0;
    }

    void consistencyY(size_t c1, size_t c2)
    {
        const LocalEdgeVector<1> S11left(
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

        tmpX.block<1, 6>(c1, 0) -= 1. / mesh.hx * avgS11 * BiG23_1;
        tmpX.block<1, 6>(c2, 0) += 1. / mesh.hx * avgS11 * BiG23_3;
        tmpY.block<1, 6>(c1, 0) -= 1. / mesh.hx * avgS12 * BiG23_1;
        tmpY.block<1, 6>(c2, 0) += 1. / mesh.hx * avgS12 * BiG23_3;
    }

    void consistencyX(size_t c1, size_t c2)
    {
        const LocalEdgeVector<1> S12bot(
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

        tmpX.block<1, 6>(c1, 0) -= 1. / mesh.hy * avgS12 * BiG23_2;
        tmpX.block<1, 6>(c2, 0) += 1. / mesh.hy * avgS12 * BiG23_0;
        tmpY.block<1, 6>(c1, 0) -= 1. / mesh.hy * avgS22 * BiG23_2;
        tmpY.block<1, 6>(c2, 0) += 1. / mesh.hy * avgS22 * BiG23_0;
    }

    template <int DGdegree>
    void boundaryDirichletLeft(const size_t c2)
    { //x=0, y=t
        const LocalEdgeVector<2> rightX(vx(c2, 0) - 0.5 * vx(c2, 1) + 1. / 6. * vx(c2, 3),
            vx(c2, 2) - 0.5 * vx(c2, 5),
            vx(c2, 4));
        const LocalEdgeVector<2> rightY(vy(c2, 0) - 0.5 * vy(c2, 1) + 1. / 6. * vy(c2, 3),
            vy(c2, 2) - 0.5 * vy(c2, 5),
            vy(c2, 4));

        tmpX.block<1, 6>(c2, 0) -= 1. / mesh.hx / mesh.hx * rightX * BiGe23 * BiG23_3;
        tmpY.block<1, 6>(c2, 0) -= 1. / mesh.hx / mesh.hx * rightY * BiGe23 * BiG23_3;
    }

    template <int DGdegree>
    void boundaryDirichletRight(const size_t c1)
    { //x=1, y=t
        const LocalEdgeVector<2> leftX(vx(c1, 0) + 0.5 * vx(c1, 1) + 1. / 6. * vx(c1, 3),
            vx(c1, 2) + 0.5 * vx(c1, 5),
            vx(c1, 4));
        const LocalEdgeVector<2> leftY(vy(c1, 0) + 0.5 * vy(c1, 1) + 1. / 6. * vy(c1, 3),
            vy(c1, 2) + 0.5 * vy(c1, 5),
            vy(c1, 4));

        tmpX.block<1, 6>(c1, 0) -= 1. / mesh.hx / mesh.hx * leftX * BiGe23 * BiG23_1;
        tmpY.block<1, 6>(c1, 0) -= 1. / mesh.hx / mesh.hx * leftY * BiGe23 * BiG23_1;
    }

    template <int DGdegree>
    void boundaryDirichletTop(const size_t c1)
    {
        const LocalEdgeVector<2> topX(
            vx(c1, 0) + 0.5 * vx(c1, 2) + 1. / 6. * vx(c1, 4),
            vx(c1, 1) + 0.5 * vx(c1, 5),
            vx(c1, 3));

        const LocalEdgeVector<2> topY(
            vy(c1, 0) + 0.5 * vy(c1, 2) + 1. / 6. * vy(c1, 4),
            vy(c1, 1) + 0.5 * vy(c1, 5),
            vy(c1, 3));

        tmpX.block<1, 6>(c1, 0) -= 1. / mesh.hy / mesh.hy * topX * BiGe23 * BiG23_2;
        tmpY.block<1, 6>(c1, 0) -= 1. / mesh.hy / mesh.hy * topY * BiGe23 * BiG23_2;
    }

    template <int DGdegree>
    void boundaryDirichletBottom(const size_t c2)
    {
        const LocalEdgeVector<2> bottomX(
            vx(c2, 0) - 0.5 * vx(c2, 2) + 1. / 6. * vx(c2, 4),
            vx(c2, 1) - 0.5 * vx(c2, 5), vx(c2, 3));

        const LocalEdgeVector<2> bottomY(
            vy(c2, 0) - 0.5 * vy(c2, 2) + 1. / 6. * vy(c2, 4),
            vy(c2, 1) - 0.5 * vy(c2, 5), vy(c2, 3));

        tmpX.block<1, 6>(c2, 0) -= 1. / mesh.hy / mesh.hy * bottomX * BiGe23 * BiG23_0;
        tmpY.block<1, 6>(c2, 0) -= 1. / mesh.hy / mesh.hy * bottomY * BiGe23 * BiG23_0;
    }

    /**!
   * controls the flow of the dynamical core
   */
    void momentumJumps();
    void momentumDirichletBoundary();

    void momentumConsistency();
    void momentumSymmetry();

    void advectionStep();
    void momentumSubsteps();
    void step();
};

}

/*----------------------------   dynamics.h     ---------------------------*/
/* end of #ifndef __dynamics_H */
#endif
/*----------------------------   dynamics.h     ---------------------------*/
