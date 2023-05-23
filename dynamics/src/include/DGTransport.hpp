/*!
 * @file    DGTransport.hpp
 * @date    Dec 5, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGTRANSPORT_HPP
#define __DGTRANSPORT_HPP

#include "cgVector.hpp"
#include "dgVector.hpp"
#include "ParametricMap.hpp"


namespace Nextsim {

template <int DG>
void DGTransportOperator(const ParametricMesh& smesh, const double dt, const DGVector<DG>& vx,
    const DGVector<DG>& vy, const EdgeVector<EDGEDOFS(DG)>& evx,
    const EdgeVector<EDGEDOFS(DG)>& evy, const DGVector<DG>& phi, DGVector<DG>& phiup);

/*!
 * Main class to manage the transport scheme on the parametric ParametricMesh
 *
 * template parameter DG, EDGEDOFS(DG) are number of local unknowns, that is
 * (1,2,3,6) on the cell and (1,2,3) on the edge
 */
template <int DG>
class DGTransport {
protected:
 
    //! spatial mesh.
    const ParametricMesh& smesh;

  //! Precomputed stencil-like matrices for efficient numerical quadrature
  ParametricTransportMap<DG> parammap;


    //! reference to the current velocity
    DGVector<DG> velx, vely;

    //! Specifies the time stepping scheme [rk1, rk2, rk3]
    std::string timesteppingscheme;

    //! normal velocity in edges parallel to X- and Y-axis
    EdgeVector<EDGEDOFS(DG)> normalvel_X, normalvel_Y;

    //! temporary vectors for time stepping
    DGVector<DG> tmp1, tmp2, tmp3;

    //! Internal functions

    /*!
     * Performs one time step transporting phi with the Fwd-Euler Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk1(const double dt, DGVector<DG>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk2(const double dt, DGVector<DG>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk3(const double dt, DGVector<DG>& phi);

public:
  DGTransport(const ParametricMesh& mesh)
      : smesh(mesh),
	parammap(mesh)	  
        , timesteppingscheme("rk2")
    {
        if (!(smesh.nelements > 0)) {
            std::cerr << "DGTransport: The mesh must already be initialized!" << std::endl;
            abort();
        }
        // Resize DG vectors for storing the velocity
        velx.resize_by_mesh(smesh);
        vely.resize_by_mesh(smesh);

        // resize tmp-vectors for time stepping
        tmp1.resize_by_mesh(smesh);
        tmp2.resize_by_mesh(smesh);
        tmp3.resize_by_mesh(smesh);

        // resize vectors to store the normal-velocity on the edges
        normalvel_Y.resize_by_mesh(smesh, EdgeType::Y);
        normalvel_X.resize_by_mesh(smesh, EdgeType::X);

	// initialize the mapping and set up required matrices
	parammap.InitializeAdvectionCellTerms();
	parammap.InitializeInverseDGMassMatrix();
    }

    // Access members

    const DGVector<DG>& GetVx() const
    {
        return velx;
    }
    const DGVector<DG>& GetVy() const
    {
        return vely;
    }
    DGVector<DG>& GetVx()
    {
        return velx;
    }
    DGVector<DG>& GetVy()
    {
        return vely;
    }

    // High level functions
    void settimesteppingscheme(const std::string tss)
    {
        timesteppingscheme = tss;
        assert((tss == "rk1") || (tss == "rk2") || (tss == "rk3"));
    }

    /*!
     * Sets the normal-velocity vector on the edges
     * The normal velocity is scaled with the length of the edge,
     * this already serves as the integraiton weight
     */
    void reinitnormalvelocity();

    /*!
     * Prepares the advection step:
     * - interpolates CG velocity to DG
     * - initializes normal velocity on the edges
     */
    template <int CG>
    void prepareAdvection(const CGVector<CG>& cg_vx, const CGVector<CG>& cg_vy);

    /*!
     * Performs one time step transporting phi
     *
     * @params phi is the vector of values to be transported
     */
    void step(const double dt, DGVector<DG>& phi);


private:
  /*!
   * Several internal functions 
   */

  //! computes all integrals on the elements and the edges
  void DGTransportOperator(const ParametricMesh& smesh, const double dt,
							  const DGVector<DG>& vx,
							  const DGVector<DG>& vy,
							  const EdgeVector<EDGEDOFS(DG)>& normalvel_X,
							  const EdgeVector<EDGEDOFS(DG)>& normalvel_Y,
							  const DGVector<DG>& phi, DGVector<DG>& phiup);

  void edge_term_X(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup, const DGVector<DG>& phi, 
		   const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie);
  void edge_term_Y(const ParametricMesh& smesh, const double dt, DGVector<DG>& phiup, const DGVector<DG>& phi, 
		   const EdgeVector<EDGEDOFS(DG)>& normalvel_Y, const size_t c1, const size_t c2, const size_t ie);
    
  void cell_term(const ParametricMesh& smesh, double dt,
		 DGVector<DG>& phiup, const DGVector<DG>& phi,
		 const DGVector<DG>& vx,
		 const DGVector<DG>& vy, const size_t ic);
  

};

} /* namespace Nextsim */

#undef EDGEDOFS

#endif /* __DGTRANSPORT_HPP */
