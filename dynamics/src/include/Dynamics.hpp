/*!
 * @file Dynamics.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DYNAMICS_HPP
#define __DYNAMICS_HPP

#include "Mesh.hpp"
#include "dgTransport.hpp"
#include "dgVector.hpp"

namespace Nextsim {
// pre-computed matrices for assembling dG-transport
// the Gauss rule for dG(n) is set to n+1 points
// this might not be enough?
#include "codeGenerationDGinGauss.hpp"

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
    CellVector<6> vx, vy; //!< velocity fields
    CellVector<3> S11, S12, S22; //!< entries of (symmetric) stress tensor
    CellVector<3> E11, E12, E22; //!< entries of (symmetric) strain stress tensor

    CellVector<6> A, H; //!< ice height and ice concentration
    CellVector<1> D; //!< ice damage. ?? Really dG(0) ??

    CellVector<6> oceanX, oceanY; //!< ocean forcing. ?? Higher order??
    CellVector<1> atmX, atmY; //!< ocean forcing. ?? Higher order??

    //! temporary vectors for time stepping
    CellVector<6> tmpX, tmpY;

    /*!
     * Subclasses for managing DG transport and the momentum problem
     */
    DGTransport<6, 3> dgtransport;

    // public:
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
    CellVector<6>& GetTMPX() { return tmpX; }
    CellVector<6>& GetTMPY() { return tmpY; }

    CellVector<6>& GetVX() { return vx; }
    CellVector<6>& GetVY() { return vy; }
    CellVector<3>& GetS11() { return S11; }
    CellVector<3>& GetS12() { return S12; }
    CellVector<3>& GetS22() { return S22; }
    CellVector<3>& GetE11() { return E11; }
    CellVector<3>& GetE12() { return E12; }
    CellVector<3>& GetE22() { return E22; }
    CellVector<1>& GetD() { return D; }
    CellVector<6>& GetH() { return H; }
    CellVector<6>& GetA() { return A; }
    // TODO: First order need to implement stabilisations etc.
    // CellVector<1>& GetVX() { return vx; }
    // CellVector<1>& GetVY() { return vy; }
    // CellVector<0>& GetS11() { return S11; }
    // CellVector<0>& GetS12() { return S12; }
    // CellVector<0>& GetS21() { return S21; }
    // CellVector<0>& GetS22() { return S22; }
    // CellVector<0>& GetD() { return D; }
    // CellVector<1>& GetH() { return H; }
    // CellVector<1>& GetA() { return A; }

    CellVector<1>& GetAtmX() { return atmX; }
    CellVector<1>& GetAtmY() { return atmY; }
    CellVector<6>& GetOceanX() { return oceanX; }
    CellVector<6>& GetOceanY() { return oceanY; }

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
    void addStressTensorBoundary(
        double scaleSigma); //!< same on boundary (from one side)
    //! consistency edge terms coming from integration by parts
    void stressTensorEdgesY(double scaleSigma, size_t c1, size_t c2)
    {
        const LocalEdgeVector<2> S11left( // left/right refers to the element
            S11(c1, 0) + 0.5 * S11(c1, 1), S11(c1, 2));
        const LocalEdgeVector<2> S11right(S11(c2, 0) - 0.5 * S11(c2, 1), S11(c2, 2));
        const LocalEdgeVector<2> S12left(S12(c1, 0) + 0.5 * S12(c1, 1), S12(c1, 2));
        const LocalEdgeVector<2> S12right(S12(c2, 0) - 0.5 * S12(c2, 1), S12(c2, 2));

        const LocalEdgeVector<3> avgS11 = 0.5 * (S11left + S11right) * PSIe<2,3>;
        const LocalEdgeVector<3> avgS12 = 0.5 * (S12left + S12right) * PSIe<2,3>;

        // {{S}} * [[Phi]]
        tmpX.block<1, 6>(c1, 0) -= scaleSigma / mesh.hx * avgS11 * PSI63_1;
        tmpX.block<1, 6>(c2, 0) += scaleSigma / mesh.hx * avgS11 * PSI63_3;
        tmpY.block<1, 6>(c1, 0) -= scaleSigma / mesh.hx * avgS12 * PSI63_1;
        tmpY.block<1, 6>(c2, 0) += scaleSigma / mesh.hx * avgS12 * PSI63_3;
    }
    void stressTensorEdgesX(double scaleSigma, size_t c1, size_t c2)
    {
        const LocalEdgeVector<2> S12bot( // bot/top refers to the element
            S12(c1, 0) + 0.5 * S12(c1, 2), S12(c1, 1));
        const LocalEdgeVector<2> S12top(S12(c2, 0) - 0.5 * S12(c2, 2), S12(c2, 1));
        const LocalEdgeVector<2> S22bot(S22(c1, 0) + 0.5 * S22(c1, 2), S22(c1, 1));
        const LocalEdgeVector<2> S22top(S22(c2, 0) - 0.5 * S22(c2, 2), S22(c2, 1));

        const LocalEdgeVector<3> avgS12 = 0.5 * (S12top + S12bot) * PSIe<2,3>;
        const LocalEdgeVector<3> avgS22 = 0.5 * (S22top + S22bot) * PSIe<2,3>;

        tmpX.block<1, 6>(c1, 0) -= scaleSigma / mesh.hy * avgS12 * PSI63_2;
        tmpX.block<1, 6>(c2, 0) += scaleSigma / mesh.hy * avgS12 * PSI63_0;
        tmpY.block<1, 6>(c1, 0) -= scaleSigma / mesh.hy * avgS22 * PSI63_2;
        tmpY.block<1, 6>(c2, 0) += scaleSigma / mesh.hy * avgS22 * PSI63_0;
    }

    //! consistency terms coming from integration by parts on boundaries
    void stressTensorBoundaryLeft(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<2> S11left( // left edge of the element
            S11(c, 0) - 0.5 * S11(c, 1), S11(c, 2));
        const LocalEdgeVector<2> S12left(S12(c, 0) - 0.5 * S12(c, 1), S12(c, 2));

        tmpX.block<1, 6>(c, 0) += scaleSigma / mesh.hx * S11left * PSIe<2,3> * PSI63_3;
        tmpY.block<1, 6>(c, 0) += scaleSigma / mesh.hx * S12left * PSIe<2,3> * PSI63_3;
    }
    void stressTensorBoundaryRight(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<2> S11right(S11(c, 0) + 0.5 * S11(c, 1), S11(c, 2));
        const LocalEdgeVector<2> S12right(S12(c, 0) + 0.5 * S12(c, 1), S12(c, 2));

        tmpX.block<1, 6>(c, 0) -= scaleSigma / mesh.hx * S11right * PSIe<2,3> * PSI63_1;
        tmpY.block<1, 6>(c, 0) -= scaleSigma / mesh.hx * S12right * PSIe<2,3> * PSI63_1;
    }
    void stressTensorBoundaryUpper(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<2> S12top(S12(c, 0) + 0.5 * S12(c, 2), S12(c, 1));
        const LocalEdgeVector<2> S22top(S22(c, 0) + 0.5 * S22(c, 2), S22(c, 1));

        tmpX.block<1, 6>(c, 0) -= scaleSigma / mesh.hy * S12top * PSIe<2,3> * PSI63_2;
        tmpY.block<1, 6>(c, 0) -= scaleSigma / mesh.hy * S22top * PSIe<2,3> * PSI63_2;
    }
    void stressTensorBoundaryLower(double scaleSigma, size_t c)
    {
        const LocalEdgeVector<2> S12bot(S12(c, 0) - 0.5 * S12(c, 2), S12(c, 1));
        const LocalEdgeVector<2> S22bot(S22(c, 0) - 0.5 * S22(c, 2), S22(c, 1));

        tmpX.block<1, 6>(c, 0) += scaleSigma / mesh.hy * S12bot * PSIe<2,3> * PSI63_0;
        tmpY.block<1, 6>(c, 0) += scaleSigma / mesh.hy * S22bot * PSIe<2,3> * PSI63_0;
    }

    //! enforce continuity of the velocity by Nitsche
public:
    void velocityContinuity(
        double gamma); //!< gamma/h < [v], [Phi] > to enforce weak continuity of velocities
private:
    void velocityContinuityY(size_t c1, size_t c2, double gamma)
    {
        const LocalEdgeVector<3> leftX(vx(c1, 0) + 0.5 * vx(c1, 1) + 1. / 6. * vx(c1, 3),
            vx(c1, 2) + 0.5 * vx(c1, 5), vx(c1, 4));
        const LocalEdgeVector<3> leftY(vy(c1, 0) + 0.5 * vy(c1, 1) + 1. / 6. * vy(c1, 3),
            vy(c1, 2) + 0.5 * vy(c1, 5), vy(c1, 4));

        const LocalEdgeVector<3> rightX(vx(c2, 0) - 0.5 * vx(c2, 1) + 1. / 6. * vx(c2, 3),
            vx(c2, 2) - 0.5 * vx(c2, 5), vx(c2, 4));
        const LocalEdgeVector<3> rightY(vy(c2, 0) - 0.5 * vy(c2, 1) + 1. / 6. * vy(c2, 3),
            vy(c2, 2) - 0.5 * vy(c2, 5), vy(c2, 4));

        const LocalEdgeVector<3> jumpX = (rightX - leftX) * PSIe<3,3>;
        tmpX.block<1, 6>(c1, 0) += gamma / mesh.hx / mesh.hx * jumpX * PSI63_1;
        tmpX.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * jumpX * PSI63_3;

        const LocalEdgeVector<3> jumpY = (rightY - leftY) * PSIe<3,3>;

        tmpY.block<1, 6>(c1, 0) += gamma / mesh.hx / mesh.hx * jumpY * PSI63_1;
        tmpY.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * jumpY * PSI63_3;
    }
    void velocityContinuityX(size_t c1, size_t c2, double gamma)
    {
        const LocalEdgeVector<3> topX(vx(c1, 0) + 0.5 * vx(c1, 2) + 1. / 6. * vx(c1, 4),
            vx(c1, 1) + 0.5 * vx(c1, 5), vx(c1, 3));
        const LocalEdgeVector<3> bottomX(vx(c2, 0) - 0.5 * vx(c2, 2) + 1. / 6. * vx(c2, 4),
            vx(c2, 1) - 0.5 * vx(c2, 5), vx(c2, 3));

        const LocalEdgeVector<3> topY(vy(c1, 0) + 0.5 * vy(c1, 2) + 1. / 6. * vy(c1, 4),
            vy(c1, 1) + 0.5 * vy(c1, 5), vy(c1, 3));
        const LocalEdgeVector<3> bottomY(vy(c2, 0) - 0.5 * vy(c2, 2) + 1. / 6. * vy(c2, 4),
            vy(c2, 1) - 0.5 * vy(c2, 5), vy(c2, 3));

        const LocalEdgeVector<3> jumpX = (topX - bottomX) * PSIe<3,3>;
        tmpX.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * jumpX * PSI63_2;
        tmpX.block<1, 6>(c2, 0) += gamma / mesh.hy / mesh.hy * jumpX * PSI63_0;

        const LocalEdgeVector<3> jumpY = (topY - bottomY) * PSIe<3,3>;
        tmpY.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * jumpY * PSI63_2;
        tmpY.block<1, 6>(c2, 0) += gamma / mesh.hy / mesh.hy * jumpY * PSI63_0;
    }

    //! weak Dirichlet boundary with Nitsche method
public:
    void velocityDirichletBoundary(double gamma);

private:
    void velocityDirichletBoundaryLeft(const size_t c2, double gamma)
    { // x=0, y=t
        const LocalEdgeVector<3> rightX(vx(c2, 0) - 0.5 * vx(c2, 1) + 1. / 6. * vx(c2, 3),
            vx(c2, 2) - 0.5 * vx(c2, 5), vx(c2, 4));
        const LocalEdgeVector<3> rightY(vy(c2, 0) - 0.5 * vy(c2, 1) + 1. / 6. * vy(c2, 3),
            vy(c2, 2) - 0.5 * vy(c2, 5), vy(c2, 4));

        tmpX.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * rightX * PSIe<3,3> * PSI63_3;
        tmpY.block<1, 6>(c2, 0) -= gamma / mesh.hx / mesh.hx * rightY * PSIe<3,3> * PSI63_3;
    }

    void velocityDirichletBoundaryRight(const size_t c1, double gamma)
    { // x=1, y=t
        const LocalEdgeVector<3> leftX(vx(c1, 0) + 0.5 * vx(c1, 1) + 1. / 6. * vx(c1, 3),
            vx(c1, 2) + 0.5 * vx(c1, 5), vx(c1, 4));
        const LocalEdgeVector<3> leftY(vy(c1, 0) + 0.5 * vy(c1, 1) + 1. / 6. * vy(c1, 3),
            vy(c1, 2) + 0.5 * vy(c1, 5), vy(c1, 4));

        tmpX.block<1, 6>(c1, 0) -= gamma / mesh.hx / mesh.hx * leftX * PSIe<3,3> * PSI63_1;
        tmpY.block<1, 6>(c1, 0) -= gamma / mesh.hx / mesh.hx * leftY * PSIe<3,3> * PSI63_1;
    }

    void velocityDirichletBoundaryTop(const size_t c1, double gamma)
    {
        const LocalEdgeVector<3> topX(vx(c1, 0) + 0.5 * vx(c1, 2) + 1. / 6. * vx(c1, 4),
            vx(c1, 1) + 0.5 * vx(c1, 5), vx(c1, 3));

        const LocalEdgeVector<3> topY(vy(c1, 0) + 0.5 * vy(c1, 2) + 1. / 6. * vy(c1, 4),
            vy(c1, 1) + 0.5 * vy(c1, 5), vy(c1, 3));

        tmpX.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * topX * PSIe<3,3> * PSI63_2;
        tmpY.block<1, 6>(c1, 0) -= gamma / mesh.hy / mesh.hy * topY * PSIe<3,3> * PSI63_2;
    }

    void velocityDirichletBoundaryBottom(const size_t c2, double gamma)
    {
        const LocalEdgeVector<3> bottomX(vx(c2, 0) - 0.5 * vx(c2, 2) + 1. / 6. * vx(c2, 4),
            vx(c2, 1) - 0.5 * vx(c2, 5), vx(c2, 3));

        const LocalEdgeVector<3> bottomY(vy(c2, 0) - 0.5 * vy(c2, 2) + 1. / 6. * vy(c2, 4),
            vy(c2, 1) - 0.5 * vy(c2, 5), vy(c2, 3));

        tmpX.block<1, 6>(c2, 0) -= gamma / mesh.hy / mesh.hy * bottomX * PSIe<3,3> * PSI63_0;
        tmpY.block<1, 6>(c2, 0) -= gamma / mesh.hy / mesh.hy * bottomY * PSIe<3,3> * PSI63_0;
    }
};

} /* namespace Nextsim */

#endif /* __DYNAMICS_HPP */
