/*!
 * @file    ParametricTransport.hpp
 * @date    July 10, 2022
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __PARAMETRICTRANSPORT_HPP
#define __PARAMETRICTRANSPORT_HPP

#include "dgVector.hpp"

namespace Nextsim {

template <int DGcell, int DGedge>
void parametricTransportOperator(const SasipMesh& smesh, const double dt, const CellVector<DGcell>& vx,
    const CellVector<DGcell>& vy, const EdgeVector<DGedge>& evx,
    const EdgeVector<DGedge>& evy, const CellVector<DGcell>& phi, CellVector<DGcell>& phiup);

/*!
 * Main class to manage the transport scheme on the parametric SasipMesh
 *
 * template parameter DGcell, DGedge are number of local unknowns, that is
 * (1,2,3,6) on the cell and (1,2,3) on the edge
 */
template <int DGcell, int DGedge>
class ParametricTransport {
protected:
    //! spatial mesh.
    const SasipMesh& smesh;

    //! reference to the current velocity
  CellVector<DGcell> velx, vely;

    //! Specifies the time stepping scheme [rk1, rk2, rk3]
    std::string timesteppingscheme;

    //! normal velocity in edges parallel to X- and Y-axis
    EdgeVector<DGedge> normalvel_X, normalvel_Y;

    //! temporary vectors for time stepping
    CellVector<DGcell> tmp1, tmp2, tmp3;

public:
  ParametricTransport(const SasipMesh& mesh)
        : smesh(mesh)
        , timesteppingscheme("rk2")
    {
      if (! (smesh.nelements>0) )
	{
	  std::cerr << "ParametricTransport: The mesh must already be initialized!" << std::endl;
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
    }


  // Access members
  
  const CellVector<DGcell>& GetVx() const
  {return velx;}
  const CellVector<DGcell>& GetVy() const
  {return vely;}
  CellVector<DGcell>& GetVx() 
  {return velx;}
  CellVector<DGcell>& GetVy() 
  {return vely;}

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
     * Performs one time step transporting phi with the Fwd-Euler Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk1(const double dt, CellVector<DGcell>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk2(const double dt, CellVector<DGcell>& phi);

    /*!
     * Performs one time step transporting phi with the 2nd Order Heun Scheme
     *
     * @params phi is the vector of values to be transported
     */
    void step_rk3(const double dt, CellVector<DGcell>& phi);

    /*!
     * Performs one time step transporting phi
     *
     * @params phi is the vector of values to be transported
     */
    void step(const double dt, CellVector<DGcell>& phi);
};

} /* namespace Nextsim */

#endif /* __PARAMETRICTRANSPORT_HPP */
