#include "ParametricMap.hpp"
#include "ParametricTools.hpp"
#include "VectorManipulations.hpp"


namespace Nextsim
{
  template<int DG>
  void ParametricTransportMap<DG>::InitializeAdvectionCellTerms()
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
	    const Eigen::Matrix<Nextsim::FloatType, 1, GP(DG)*GP(DG)> cos_lat = (ParametricTools::getGaussPointsInElement<GP(DG)>(smesh, eid).row(1).array()).cos();
	    AdvectionCellTermY[eid] = PSIy<DG, GP(DG)>.array().rowwise() * (dxT.row(0).array() * cos_lat.array()) - PSIx<DG, GP(DG)>.array().rowwise() * (dyT.row(0).array() * cos_lat.array());
	  }
	else if (type == CARTESIAN)
	  AdvectionCellTermY[eid] = PSIy<DG, GP(DG)>.array().rowwise() * dxT.row(0).array() - PSIx<DG, GP(DG)>.array().rowwise() * dyT.row(0).array();
	else abort();
      }
  }



  template<int DG>
  void ParametricTransportMap<DG>::InitializeInverseDGMassMatrix()
    {
    // for advection
      InverseDGMassMatrix.clear();
      InverseDGMassMatrix.resize(smesh.nelements);


      if (type == SPHERICAL)
	{
#pragma omp parallel for
	  for (size_t eid = 0; eid<smesh.nelements;++eid)
	    InverseDGMassMatrix[eid] = SphericalTools::massMatrix<DG>(smesh, eid).inverse() / Nextsim::EarthRadius;
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


  //////////////////////////////////////////////////
  ////////////////////////////////////////////////// Momentum
  //////////////////////////////////////////////////

  
  //!
  template<int CG>
  void ParametricMomentumMap<CG>::InitializeLumpedCGMassMatrix()
  {
    lumpedcgmass.resize_by_mesh(smesh);

    for (size_t i = 0; i < smesh.nnodes; ++i)
      lumpedcgmass(i, 0) = 0;

#define CGGP(CG) ( (CG==1?3:4) )

    
    for (size_t p=0;p<2;++p) // for parallelization
      {
#pragma omp parallel for
	for (size_t iy=p;iy<smesh.ny;iy+=2)
	  for (size_t ix=0;ix<smesh.nx;++ix)
	    {
	      size_t eid = smesh.ny*iy+ix;


	      Eigen::Vector<Nextsim::FloatType, (CG==1?4:9) > Meid;

	      if (type == CARTESIAN)
		{
		  const Eigen::Matrix<Nextsim::FloatType, 1, CGGP(CG) * CGGP(CG) > J
		    = ParametricTools::J<CGGP(CG)>(smesh, eid).array() * GAUSSWEIGHTS<CGGP(CG) >.array();
		  
		  Meid = PHI<CG,CGGP(CG)> * J.transpose();
		}
	      else if (type == SPHERICAL)
		{
		  const Eigen::Matrix<Nextsim::FloatType, 1, CGGP(CG)*CGGP(CG)> cos_lat = (ParametricTools::getGaussPointsInElement<CGGP(CG)>(smesh, eid).row(1).array()).cos();

		  const Eigen::Matrix<Nextsim::FloatType, 1, CGGP(CG) * CGGP(CG) > J
		    = ParametricTools::J<CGGP(CG)>(smesh, eid).array() * GAUSSWEIGHTS<CGGP(CG) >.array() * cos_lat.array();
		  
		  Meid = PHI<CG,CGGP(CG)> * J.transpose();
		}
	      else abort();
	      
	      // index of first dof in element
	      const size_t sy = CG * smesh.nx + 1;
	      const size_t n0 = CG * iy * sy + CG * ix;
	      
	      if (CG==1)
		{
		  lumpedcgmass(n0, 0)                += Meid(0);
		  lumpedcgmass(n0 + 1, 0)            += Meid(1);
		  
		  lumpedcgmass(n0     + sy, 0)       += Meid(2);
		  lumpedcgmass(n0 + 1 + sy, 0)       += Meid(3);
		}
	      else if (CG==2)
		{
		  lumpedcgmass(n0, 0)                += Meid(0);
		  lumpedcgmass(n0 + 1, 0)            += Meid(1);
		  lumpedcgmass(n0 + 2, 0)            += Meid(2);
		  
		  lumpedcgmass(n0     + sy, 0)       += Meid(3);
		  lumpedcgmass(n0 + 1 + sy, 0)       += Meid(4);
		  lumpedcgmass(n0 + 2 + sy, 0)       += Meid(5);
		  
		  lumpedcgmass(n0     + 2 * sy, 0)   += Meid(6);
		  lumpedcgmass(n0 + 1 + 2 * sy, 0)   += Meid(7);
		  lumpedcgmass(n0 + 2 + 2 * sy, 0)   += Meid(8);
		}
	      else
		abort();
	    }
      }
    VectorManipulations::CGAddPeriodic(smesh, lumpedcgmass);
  }



    //!
  template<int CG>
  void ParametricMomentumMap<CG>::InitializeDivSMatrices()
  {
    divS1.resize(smesh.nelements);
    divS2.resize(smesh.nelements);

    // parallel loop over all elements for computing entries
#pragma omp parallel for
    for (size_t eid = 0; eid < smesh.nelements; ++eid) {
      
      //               [  dyT2  -dxT2 ]
      // A = JF^{-T} = [              ]
      //               [ -dyT1   dxT1 ]
      //
      const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(CG2DGSTRESS(CG))> dxT = (ParametricTools::dxT<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array()).matrix();
      const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(CG2DGSTRESS(CG))> dyT = (ParametricTools::dyT<GAUSSPOINTS1D(CG2DGSTRESS(CG))>(smesh, eid).array().rowwise() * GAUSSWEIGHTS<GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array()).matrix();
      
      // the transformed gradient of the CG basis function in the gauss points
      const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), GAUSSPOINTS(CG2DGSTRESS(CG))> dx_cg2 = PHIx<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * dyT.row(1).array() - PHIy<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * dxT.row(1).array();
      
      const Eigen::Matrix<Nextsim::FloatType, (CG == 2 ? 9 : 4), GAUSSPOINTS(CG2DGSTRESS(CG))> dy_cg2 = PHIy<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * dxT.row(0).array() - PHIx<CG, GAUSSPOINTS1D(CG2DGSTRESS(CG))>.array().rowwise() * dyT.row(0).array();
      
      // PSI83 is the DG-Basis-function in the guass-point
      // PSI83_{iq} = PSI_i(q)   [ should be called DG_in_GAUSS ]

      
      divS1[eid] = dx_cg2 * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose();
      divS2[eid] = dy_cg2 * PSI<CG2DGSTRESS(CG), GAUSSPOINTS1D(CG2DGSTRESS(CG))>.transpose();
    }

  }



  
  template class ParametricTransportMap<1>;
  template class ParametricTransportMap<3>;
  template class ParametricTransportMap<6>;

    
  template class ParametricMomentumMap<1>;
    template class ParametricMomentumMap<2>;
  

}


