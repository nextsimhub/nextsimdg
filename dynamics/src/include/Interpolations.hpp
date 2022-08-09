/*!
 * @file Interpolations.hpp
 * @date Aug 9 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __INTERPOLATIONS_HPP
#define __INTERPOLATIONS_HPP

#include "dgVector.hpp"
#include "cgVector.hpp"


/*!
 * This file includes functions to interpolate between different vectors
 * that are used in the dynamics core, advection and momentum
 * 
 * - Function2CG: Interpolates a function to a CG vector
 * - Function2DG: Projects a function to a DG vector
 * - CG2DG:       Projects a CG vector to a DG vector
 * - DG2GG:       Interpolates a DG Vector to a CG vector
 *
 * Interpolations just pick the node-wise values, where an interpolation from
 * DG to CG will average the values on the edges and the vertices
 *
 * Projections are performed as L2-projectio such that the local
 * DG mass matrix must in inverted. 
 */

namespace Nextsim {


  namespace Interpolations {

    //! interface class to provide an analytical function
    class Function
    {
    public:
      virtual double operator()(double x, double y) const=0;
    };


    template<int CG>
    void Function2CG(const SasipMesh& smesh, CGVector<CG>& dest, const Function& src);
    template<int DG>
    void Function2DG(const SasipMesh& smesh, CellVector<DG>& dest, const Function& src);
    template<int CG, int DG>
    void CG2DG(const SasipMesh& smesh, CellVector<DG>& dest, const CGVector<CG> &src);
    template<int CG, int DG>
    void DG2CG(const SasipMesh& smesh, CGVector<CG>& dest, const CellVector<DG> &src);
    
  } /* namespace Interpolation */

} /* namespace Nextsim */

#endif /* __INTERPOLATIONS_HPP */
