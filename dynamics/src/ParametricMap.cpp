#include "ParametricMap.hpp"
#include "ParametricTools.hpp"

namespace Nextsim
{
  template<int DG>
  void ParametricMap<DG>::InitializeAdvectionCellTerms()
    {
    // for advection
    AdvectionCellTermX.clear();
    AdvectionCellTermY.clear();
    
    AdvectionCellTermX.resize(smesh.nelements);
    AdvectionCellTermY.resize(smesh.nelements);

        //    phiup.row(eid) += dt * (parammap.AdvectionCellTermX[eid].array().rowwise() * vx_gauss.array() + parammap.AdvectionCellTermY[eid].array().rowwise() * vy_gauss.array()).matrix() * phi_gauss.transpose();
    
    // gradient of transformation
    //      [ dxT1, dyT1 ]     //            [ dyT2, -dxT2 ]
    // dT = 		     // J dT^{-T}=
    //      [ dxT2, dyT2 ]     //            [ -dyT1, dxT1 ]
    //
    // given as [dxT1, dxT2, dyT1, dyT2] ->  [dyT2, -dxT2, -dyT1, dxT1 ]

    // J dT^{-T} nabla Phi  = [dyT2 * PSIx - dxT2 * PSIy, -dyT1 * PSIx + dxT1 * PSIy]
    // PSIx, PSIy are DG x QQ - matrices
    // dxT, dyT are 2 x QQ - matrices

    // Store wq * phi(q)

    
#pragma omp parallel for
    for (size_t eid = 0; eid<smesh.nelements;++eid)
      {
	const Eigen::Matrix<Nextsim::FloatType, 2, GP(DG)* GP(DG)> dxT = ParametricTools::dxT<GP(DG)>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<GP(DG)>.array();
	const Eigen::Matrix<Nextsim::FloatType, 2, GP(DG)* GP(DG)> dyT = ParametricTools::dyT<GP(DG)>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<GP(DG)>.array();
	
	// [J dT^{-T} nabla phi]_1
	AdvectionCellTermX[eid] = PSIx<DG, GP(DG)>.array().rowwise() * dyT.row(1).array() - PSIy<DG, GP(DG)>.array().rowwise() * dxT.row(1).array();
	// [J dT^{-T} nabla phi]_2

	//! the lat-direction must be scaled with the metric term if in the spherical system
	if (type == SPHERICAL)
	  {
	    const Eigen::Matrix<Nextsim::FloatType, 1, GP(DG)*GP(DG)> cos_lat = (ParametricTools::getGaussPointsInElement<GP(DG)>(smesh, eid).row(1).array()*M_PI/180.).cos();
	    AdvectionCellTermY[eid] = PSIy<DG, GP(DG)>.array().rowwise() * (dxT.row(0).array() * cos_lat.array()) - PSIx<DG, GP(DG)>.array().rowwise() * (dyT.row(0).array() * cos_lat.array());
	  }
	else if (type == CARTESIAN)
	  AdvectionCellTermY[eid] = PSIy<DG, GP(DG)>.array().rowwise() * dxT.row(0).array() - PSIx<DG, GP(DG)>.array().rowwise() * dyT.row(0).array();
	else abort();
      }
  }



  template<int DG>
  void ParametricMap<DG>::InitializeInverseDGMassMatrix()
    {
    // for advection
      InverseDGMassMatrix.clear();
      InverseDGMassMatrix.resize(smesh.nelements);

      
      if (type == SPHERICAL)
	{
#pragma omp parallel for
	  for (size_t eid = 0; eid<smesh.nelements;++eid)
	    InverseDGMassMatrix[eid] = SphericalTools::massMatrix<DG>(smesh, eid).inverse();
	}
      else if (type == CARTESIAN)
	{
#pragma omp parallel for
	  for (size_t eid = 0; eid<smesh.nelements;++eid)
	    InverseDGMassMatrix[eid] = ParametricTools::massMatrix<DG>(smesh, eid).inverse();
	}
      else
	{
	  std::cerr << "Coordinate System " << type << " not known!" << std::endl;
	  abort();
	}
  }

  template class ParametricMap<1>;
  template class ParametricMap<3>;
  template class ParametricMap<6>;


}


