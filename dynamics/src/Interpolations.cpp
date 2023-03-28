#include "Interpolations.hpp"
#include "ParametricTools.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"

namespace Nextsim {
  namespace Interpolations {

    // ******************** Function -> CG ******************** //

    template <int CG>
    void Function2CG(const ParametricMesh& smesh, CGVector<CG>& dest, const Function& src)
    {
      assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == dest.rows());

#pragma omp parallel for
      for (size_t iy = 0; iy < CG * smesh.ny + 1; ++iy) {
	size_t ii = iy * (CG * smesh.nx + 1); // first index in row

	for (size_t ix = 0; ix < CG * smesh.nx + 1; ++ix, ++ii) {
	  const auto v = smesh.coordinate<CG>(ix, iy);
	  dest(ii) = src(v[0], v[1]);
	}
      }
    }

    // ******************** Function -> DG ******************** //



    // ******************** Function -> DG ******************** //

    template <>
    void Function2DG(const ParametricMesh& smesh, DGVector<1>& phi, const Function& initial)
    { // 2point-gauss rule
      phi.setZero();

#pragma omp parallel for
      for (size_t eid = 0; eid < smesh.nelements; ++eid) {
	if (!smesh.landmask[eid])
	  continue;
	const double mass = smesh.area(eid);
	// transform gauss points to real element

	// GAUSSPOINTS_3 is the 2 x 4 - Matrix with the reference coordinates
	// BIG33 is the 3 x 4 matrix with DG-Psi_i(q), i=0,1,2; q=0,1,...,8
	// CG_Q1_3 is the 4 x 4 matrix with the four Q1-basis functions in the GP
	// coordinates is the 4 x 2 - matrix with the coords of the 4 vertices

	// the Gauss points in the element 2 x 4 - Matrix

	const Eigen::Matrix<Nextsim::FloatType, 2, 4> gp     = ParametricTools::getGaussPointsInElement<2>(smesh, eid);
	const Eigen::Matrix<Nextsim::FloatType, 4, 1>
	  initial_in_gp(initial(gp(0, 0), gp(1, 0)),
			initial(gp(0, 1), gp(1, 1)),
			initial(gp(0, 2), gp(1, 2)),
			initial(gp(0, 3), gp(1, 3)));
	const Eigen::Matrix<Nextsim::FloatType, 1, 4> J      = ParametricTools::J<2>(smesh,eid);

	if (smesh.CoordinateSystem == SPHERICAL)
	  {
	    const Eigen::Matrix<Nextsim::FloatType, 1, 4> coslat = (gp.row(1).array()).cos();
	    const double masscosJ =  (J.array() * GAUSSWEIGHTS<2>.array() * coslat.array()).sum();
	    phi.row(eid) = 1. / masscosJ * (((coslat.array() * ParametricTools::J<2>(smesh, eid).array() * GAUSSWEIGHTS<2>.array())).matrix() * initial_in_gp);
	  }
	else if (smesh.CoordinateSystem == CARTESIAN)
	  {
	    const double masscosJ =  (J.array() * GAUSSWEIGHTS<2>.array()).sum();
	    phi.row(eid) = 1. / masscosJ * (((ParametricTools::J<2>(smesh, eid).array() * GAUSSWEIGHTS<2>.array())).matrix() * initial_in_gp);
	  }
	else abort();
      }
    }

    template <int DG>
    void Function2DG(const ParametricMesh& smesh, DGVector<DG>& phi, const Function& initial)
    { // 2point-gauss rule
      phi.setZero();

#pragma omp parallel for
      for (size_t eid = 0; eid < smesh.nelements; ++eid) {
	if (!smesh.landmask[eid])
	  continue;

	// transform gauss points to real element

	// GAUSSPOINTS_3 is the 2 x 4 - Matrix with the reference coordinates
	// BIG33 is the 3 x 4 matrix with DG-Psi_i(q), i=0,1,2; q=0,1,...,8
	// CG_Q1_3 is the 4 x 4 matrix with the four Q1-basis functions in the GP
	// coordinates is the 4 x 2 - matrix with the coords of the 4 vertices

	// the Gauss points in the element 2 x 4 - Matrix
	const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS(DG)> gp     = ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(smesh, eid);
	Eigen::Matrix<Nextsim::FloatType, GAUSSPOINTS(DG), 1> initial_in_gp;
	for (size_t i=0;i<GAUSSPOINTS1D(DG)*GAUSSPOINTS1D(DG);++i)
	  initial_in_gp(i) = initial(gp(0, i), gp(1, i));

	if (smesh.CoordinateSystem == SPHERICAL)
	  {
	    const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS(DG)> coslat = (gp.row(1).array()).cos();
	    phi.row(eid) = SphericalTools::massMatrix<DG>(smesh,eid).inverse() * ((PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, eid).array() * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array() * coslat.array())).matrix() * initial_in_gp);
	  }
	else if (smesh.CoordinateSystem == CARTESIAN)
	  {
	    phi.row(eid) = ParametricTools::massMatrix<DG>(smesh,eid).inverse() * ((PSI<DG, GAUSSPOINTS1D(DG)>.array().rowwise() * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, eid).array() * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array())).matrix() * initial_in_gp);
	  }
	else abort();
      }
    }


    template <int CG, int DG>
    void CG2DG(const ParametricMesh& smesh, DGVector<DG>& dg, const CGVector<CG>& cg)
    { 
      assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == cg.rows());
      assert(static_cast<long int>(smesh.nx * smesh.ny) == dg.rows());

      const int cgshift = CG * smesh.nx + 1; //!< Index shift for each row

      // parallelize over elements
#pragma omp parallel for
      for (size_t dgi = 0; dgi < smesh.nelements; ++dgi) {
	size_t iy = dgi / smesh.nx; //!< y-index of element
	size_t ix = dgi % smesh.nx; //!< x-index of element

	size_t cgi = CG * cgshift * iy + CG * ix; //!< lower/left Index in cg vector

	Eigen::Matrix<double, (CG == 2 ? 9 : 4), 1> cg_local; //!< the 9 local unknowns in the element
	if (CG == 1) {
	  cg_local << cg(cgi), cg(cgi + 1), cg(cgi + cgshift), cg(cgi + 1 + cgshift);
	} else {
	  cg_local << cg(cgi), cg(cgi + 1), cg(cgi + 2), cg(cgi + cgshift), cg(cgi + 1 + cgshift),
	    cg(cgi + 2 + cgshift), cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift),
	    cg(cgi + 2 + 2 * cgshift);
	}
	// solve:  (Vdg, PHI) = (Vcg, PHI) with mapping to spher. coord.
	if (smesh.CoordinateSystem == SPHERICAL)
	  dg.row(dgi) =
	    SphericalTools::massMatrix<DG>(smesh, dgi).inverse() * PSI<DG, GAUSSPOINTS1D(DG)>
	    * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, dgi).array() *
	       GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array() *
	       (ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(smesh, dgi).row(1).array()).cos() * //! metric term
	       (PHI<CG, GAUSSPOINTS1D(DG)>.transpose() * cg_local).transpose().array()).matrix().transpose();
	else if (smesh.CoordinateSystem == CARTESIAN)
	  dg.row(dgi) =
	    ParametricTools::massMatrix<DG>(smesh, dgi).inverse() * PSI<DG, GAUSSPOINTS1D(DG)>
	    * (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, dgi).array() *
	       GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array() *
	       (PHI<CG, GAUSSPOINTS1D(DG)>.transpose() * cg_local).transpose().array()).matrix().transpose();


      }
    }

    // ******************** DG -> CG  ******************** //

    // For DG2CG, the kind of coordinate system should not matter. This is no projection
    // but just an interpolation. And this is done on the level of the reference coordinate
    // system by restricting the test functions to the lagrange points.
    template <int DG>
    void DG2CGCell(const ParametricMesh& smesh, const size_t c, const size_t cx, const size_t cy,
		   CGVector<1>& cg_A, const DGVector<DG>& A)
    {
      const size_t CGDofsPerRow = smesh.nx + 1;
      const size_t cgi = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)

      const Eigen::Matrix<double, 1, 4> At = A.row(c) * PSILagrange<DG, 2>;

      cg_A(cgi) += 0.25 * At(0);
      cg_A(cgi + 1) += 0.25 * At(1);
      cg_A(cgi + CGDofsPerRow) += 0.25 * At(2);
      cg_A(cgi + CGDofsPerRow + 1) += 0.25 * At(3);
    }
    template <int DG>
    void DG2CGCell(const ParametricMesh& smesh, const size_t c, const size_t cx, const size_t cy,
		   CGVector<2>& cg_A, const DGVector<DG>& A)
    {
      const size_t CGDofsPerRow = 2 * smesh.nx + 1;
      const size_t cgi = 2 * CGDofsPerRow * cy
	+ 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)

      const Eigen::Matrix<double, 1, 9> At = A.row(c) * PSILagrange<DG, 3>;

      cg_A(cgi) += 0.25 * At(0);
      cg_A(cgi + 1) += 0.5 * At(1);
      cg_A(cgi + 2) += 0.25 * At(2);
      cg_A(cgi + CGDofsPerRow) += 0.5 * At(3);
      cg_A(cgi + CGDofsPerRow + 1) += At(4);
      cg_A(cgi + CGDofsPerRow + 2) += 0.5 * At(5);
      cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * At(6);
      cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * At(7);
      cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * At(8);
    }

    template<int CG>
    void DG2CGBoundary(const ParametricMesh& smesh, CGVector<CG>& cg_A)
    {
      const size_t CGDofsPerRow = CG * smesh.nx + 1;
      const size_t UpperLeftIndex = CG * CGDofsPerRow * smesh.ny;

#pragma omp parallel for
      for (size_t i = 0; i < CG * smesh.nx + 1; ++i) {
	cg_A(i) *= 2.0;
	cg_A(UpperLeftIndex + i) *= 2.0;
      }
#pragma omp parallel for
      for (size_t i = 0; i < CG * smesh.ny + 1; ++i) {
	cg_A(i * CGDofsPerRow) *= 2.0;
	cg_A(i * CGDofsPerRow + CG * smesh.nx) *= 2.0;
      }
    }

    template <int CG, int DG>
    void DG2CG(const ParametricMesh& smesh, CGVector<CG>& dest, const DGVector<DG>& src)
    {
      assert(src.rows() == static_cast<int>(smesh.nx * smesh.ny));
      assert(dest.rows() == static_cast<int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)));

      dest.zero();

      // parallelization by running over stripes
      for (size_t p = 0; p < 2; ++p) {
#pragma omp parallel for
	for (size_t cy = 0; cy < smesh.ny; ++cy) {
	  if (cy % 2 == p)
	    continue;

	  size_t c = cy * smesh.nx;

	  for (size_t cx = 0; cx < smesh.nx; ++cx, ++c)
	    DG2CGCell(smesh, c, cx, cy, dest, src);
	}
      }
      DG2CGBoundary(smesh, dest);
    }


    template <int DG>
    double L2ErrorFunctionDG(const ParametricMesh& smesh, const DGVector<DG>& src, const Function& initial)
    {
      double error = 0;

#pragma omp parallel for reduction(+		\
                                   : error)
      for (size_t eid = 0; eid < smesh.nelements; ++eid) {
	if (!smesh.landmask[eid])
	  continue;

	const Eigen::Matrix<Nextsim::FloatType, 2, GAUSSPOINTS1D(DG)* GAUSSPOINTS1D(DG)> gp = ParametricTools::getGaussPointsInElement<GAUSSPOINTS1D(DG)>(smesh, eid);
	const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS1D(DG)*GAUSSPOINTS1D(DG)> src_in_gauss = src.row(eid) * PSI<DG, GAUSSPOINTS1D(DG)>;

	Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS1D(DG)*GAUSSPOINTS1D(DG)> initial_in_gp;
	for (size_t i=0;i<GAUSSPOINTS1D(DG)*GAUSSPOINTS1D(DG);++i)
	  initial_in_gp(i) = initial(gp(0,i), gp(1,i));

	// Jq * wq * Psi_i(x_q) * f(x_q)
	// matrix of size 3 x 4

	if (smesh.CoordinateSystem == SPHERICAL)
	  {
	    const Eigen::Matrix<Nextsim::FloatType, 1, GAUSSPOINTS1D(DG)*GAUSSPOINTS1D(DG)> cos_lat =
	      (gp.row(1).array()).cos();
	    error += (cos_lat.array() * ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, eid).array() * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array() * (src_in_gauss.array() - initial_in_gp.array()).square()).sum();
	  }
	else if (smesh.CoordinateSystem == CARTESIAN)
	  error += (ParametricTools::J<GAUSSPOINTS1D(DG)>(smesh, eid).array() * GAUSSWEIGHTS<GAUSSPOINTS1D(DG)>.array() * (src_in_gauss.array() - initial_in_gp.array()).square()).sum();
	else abort();
      }
      return error;
    }



    template <int CG>
    double L2ErrorFunctionCG(const ParametricMesh& smesh, const CGVector<CG>& src, const Function& initial)
    {
      double error = 0;

      const int cgshift = CG * smesh.nx + 1; //!< Index shift for each row

      // parallelize over the rows
#pragma omp parallel for reduction(+		\
                                   : error)
      for (size_t row = 0; row < smesh.ny; ++row) {
	size_t eid = smesh.nx * row; //!< Index of dg vector
	size_t cgi = CG * cgshift * row; //!< Lower left index of cg vector

	for (size_t col = 0; col < smesh.nx; ++col, ++eid, cgi += CG) { // loop over all elements

	  if (smesh.landmask[eid]==0) // only on ice
	    continue;


	  // get the local CG unknowns
	  Eigen::Matrix<double, CGDOFS(CG), 1> local;
	  if (CG == 1) {
	    local << src(cgi), src(cgi + 1), src(cgi + cgshift), src(cgi + 1 + cgshift);
	  } else if (CG == 2) {
	    local << src(cgi), src(cgi + 1), src(cgi + 2), src(cgi + cgshift), src(cgi + 1 + cgshift),
	      src(cgi + 2 + cgshift), src(cgi + 2 * cgshift), src(cgi + 1 + 2 * cgshift),
	      src(cgi + 2 + 2 * cgshift);
	  } else
	    abort();


	  const Eigen::Matrix<Nextsim::FloatType, 2, 3*3> gp = ParametricTools::getGaussPointsInElement<3>(smesh, eid);
	  const Eigen::Matrix<Nextsim::FloatType, 1, 3*3> src_in_gauss = local.transpose() * PHI<CG, 3>;

	  Eigen::Matrix<Nextsim::FloatType, 1, 3*3> initial_in_gp;
	  for (size_t i=0;i<3*3;++i)
	    initial_in_gp(i) = initial(gp(0,i), gp(1,i));

	  if (smesh.CoordinateSystem == SPHERICAL)
	    {
	      const Eigen::Matrix<Nextsim::FloatType, 1, 3*3> cos_lat =
		(gp.row(1).array()).cos();
	      error += (cos_lat.array()
			* ParametricTools::J<3>(smesh, eid).array()
			* GAUSSWEIGHTS<3>.array()
			* (src_in_gauss.array() - initial_in_gp.array()).square()).sum();
	    }
	  else if (smesh.CoordinateSystem == CARTESIAN)
	    error += (ParametricTools::J<3>(smesh, eid).array()
		      * GAUSSWEIGHTS<3>.array()
		      * (src_in_gauss.array() - initial_in_gp.array()).square()).sum();
	  else abort();
	}
      }
      return error;
    }


    template void DG2CG(const ParametricMesh& smesh, CGVector<2>& dest, const DGVector<1>& src);
    template void DG2CG(const ParametricMesh& smesh, CGVector<2>& dest, const DGVector<3>& src);
    template void DG2CG(const ParametricMesh& smesh, CGVector<2>& dest, const DGVector<6>& src);
    template void DG2CG(const ParametricMesh& smesh, CGVector<1>& dest, const DGVector<1>& src);
    template void DG2CG(const ParametricMesh& smesh, CGVector<1>& dest, const DGVector<3>& src);
    template void DG2CG(const ParametricMesh& smesh, CGVector<1>& dest, const DGVector<6>& src);

    template void CG2DG(const ParametricMesh& smesh, DGVector<1>& dg, const CGVector<1>& cg);
    template void CG2DG(const ParametricMesh& smesh, DGVector<3>& dg, const CGVector<1>& cg);
    template void CG2DG(const ParametricMesh& smesh, DGVector<6>& dg, const CGVector<1>& cg);
    template void CG2DG(const ParametricMesh& smesh, DGVector<1>& dg, const CGVector<2>& cg);
    template void CG2DG(const ParametricMesh& smesh, DGVector<3>& dg, const CGVector<2>& cg);
    template void CG2DG(const ParametricMesh& smesh, DGVector<6>& dg, const CGVector<2>& cg);

    template void Function2CG(const ParametricMesh& smesh, CGVector<1>& phi, const Function& initial);
    template void Function2CG(const ParametricMesh& smesh, CGVector<2>& phi, const Function& initial);

    template void Function2DG(const ParametricMesh& smesh, DGVector<3>& phi, const Function& initial);
    template void Function2DG(const ParametricMesh& smesh, DGVector<6>& phi, const Function& initial);
    template void Function2DG(const ParametricMesh& smesh, DGVector<8>& phi, const Function& initial);

    template double L2ErrorFunctionDG(const ParametricMesh& smesh, const DGVector<1>& src, const Function& fct);
    template double L2ErrorFunctionDG(const ParametricMesh& smesh, const DGVector<3>& src, const Function& fct);
    template double L2ErrorFunctionDG(const ParametricMesh& smesh, const DGVector<6>& src, const Function& fct);
    template double L2ErrorFunctionDG(const ParametricMesh& smesh, const DGVector<8>& src, const Function& fct);

    template double L2ErrorFunctionCG(const ParametricMesh& smesh, const CGVector<1>& src, const Function& fct);
    template double L2ErrorFunctionCG(const ParametricMesh& smesh, const CGVector<2>& src, const Function& fct);

  }
}
