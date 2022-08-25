#include "Interpolations.hpp"
#include "ParametricTools.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "codeGenerationDGinGauss.hpp"

namespace Nextsim {
namespace Interpolations {

    // ******************** Function -> CG ******************** //
    template <>
    void Function2CG(const SasipMesh& smesh, CGVector<1>& dest, const Function& src)
    {
        assert(static_cast<long int>((smesh.nx + 1) * (smesh.ny + 1)) == dest.rows());

#pragma omp parallel for
        for (size_t iy = 0; iy < smesh.ny + 1; ++iy) {
            size_t ii = iy * (smesh.nx + 1);

            for (size_t ix = 0; ix < smesh.nx + 1; ++ix, ++ii) {
                const auto v = smesh.coordinate<1>(ix, iy);
                dest(ii) = src(v[0], v[1]);
            }
        }
    }

    template <>
    void Function2CG(const SasipMesh& smesh, CGVector<2>& dest, const Function& src)
    {
        assert(static_cast<long int>((2 * smesh.nx + 1) * (2 * smesh.ny + 1)) == dest.rows());

#pragma omp parallel for
        for (size_t iy = 0; iy < 2 * smesh.ny + 1; ++iy) {
            size_t ii = iy * (2 * smesh.nx + 1);

            for (size_t ix = 0; ix < 2 * smesh.nx + 1; ++ix, ++ii) {
                const auto v = smesh.coordinate<2>(ix, iy);
                dest(ii) = src(v[0], v[1]);
            }
        }
    }

    // ******************** Function -> DG ******************** //

    template <>
    void Function2DG(const SasipMesh& smesh, CellVector<1>& phi, const Function& initial)
    { // 2point-gauss rule
        phi.setZero();

#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid) {
            const double mass = smesh.area(eid);
            // transform gauss points to real element

            // GAUSSPOINTS_3 is the 2 x 4 - Matrix with the reference coordinates
            // BIG33 is the 3 x 4 matrix with DG-Psi_i(q), i=0,1,2; q=0,1,...,8
            // CG_Q1_3 is the 4 x 4 matrix with the four Q1-basis functions in the GP
            // coordinates is the 4 x 2 - matrix with the coords of the 4 vertices

            // the Gauss points in the element 2 x 4 - Matrix
            const auto gp = ParametricTools::getGaussPointsInElement2(smesh, eid);

            const Eigen::Matrix<Nextsim::FloatType, 4, 1>
                initial_in_gp(initial(gp(0, 0), gp(1, 0)),
                    initial(gp(0, 1), gp(1, 1)),
                    initial(gp(0, 2), gp(1, 2)),
                    initial(gp(0, 3), gp(1, 3)));

            // Jq * wq * Psi_i(x_q) * f(x_q)
            // matrix of size 3 x 4
	    phi.row(eid) = 1. / mass * (((ParametricTools::J<2>(smesh, eid).array() * GAUSSWEIGHTS<2>.array())).matrix() * initial_in_gp);
        }
    }

    template <>
    void Function2DG(const SasipMesh& smesh, CellVector<3>& phi, const Function& initial)
    { // 2point-gauss rule
        phi.setZero();

#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid) {
            const Eigen::Matrix<Nextsim::FloatType, 3, 3> mass
                = Nextsim::ParametricTools::massMatrix<3>(smesh, eid);
            // transform gauss points to real element

            // GAUSSPOINTS_3 is the 2 x 4 - Matrix with the reference coordinates
            // BIG33 is the 3 x 4 matrix with DG-Psi_i(q), i=0,1,2; q=0,1,...,8
            // CG_Q1_3 is the 4 x 4 matrix with the four Q1-basis functions in the GP
            // coordinates is the 4 x 2 - matrix with the coords of the 4 vertices

            // the Gauss points in the element 2 x 4 - Matrix
            const auto gp = ParametricTools::getGaussPointsInElement2(smesh, eid);

            const Eigen::Matrix<Nextsim::FloatType, 4, 1>
                initial_in_gp(initial(gp(0, 0), gp(1, 0)),
                    initial(gp(0, 1), gp(1, 1)),
                    initial(gp(0, 2), gp(1, 2)),
                    initial(gp(0, 3), gp(1, 3)));

            // Jq * wq * Psi_i(x_q) * f(x_q)
            // matrix of size 3 x 4
	    phi.row(eid) = mass.inverse() * ((PSI<3,2>.array().rowwise() * (ParametricTools::J<2>(smesh, eid).array() * GAUSSWEIGHTS<2>.array())).matrix() * initial_in_gp);
        }
    }

    template <>
    void Function2DG(const SasipMesh& smesh, CellVector<6>& phi, const Function& initial)
    { // 3point-gauss rule
        phi.setZero();

#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid) {
            const Eigen::Matrix<Nextsim::FloatType, 6, 6> mass
                = Nextsim::ParametricTools::massMatrix<6>(smesh, eid);
            // transform gauss points to real element

            // GAUSSPOINTS_3 is the 2 x 9 - Matrix with the reference coordinates
            // BIG33 is the 3 x 9 matrix with DG-Psi_i(q), i=0,1,2; q=0,1,...,8
            // CG_Q1_3 is the 4 x 9 matrix with the four Q1-basis functions in the GP
            // coordinates is the 4 x 2 - matrix with the coords of the 4 vertices

            // the Gauss points in the element 2 x 9 - Matrix
            const auto gp = ParametricTools::getGaussPointsInElement3(smesh, eid);

            const Eigen::Matrix<Nextsim::FloatType, 9, 1>
                initial_in_gp({ { initial(gp(0, 0), gp(1, 0)),
                    initial(gp(0, 1), gp(1, 1)),
                    initial(gp(0, 2), gp(1, 2)),
                    initial(gp(0, 3), gp(1, 3)),
                    initial(gp(0, 4), gp(1, 4)),
                    initial(gp(0, 5), gp(1, 5)),
                    initial(gp(0, 6), gp(1, 6)),
                    initial(gp(0, 7), gp(1, 7)),
                    initial(gp(0, 8), gp(1, 8)) } });

            // Jq * wq * Psi_i(x_q) * f(x_q)
            // matrix of size 3 x 9
	    phi.row(eid) = mass.inverse() * ((PSI<6,3>.array().rowwise() * (ParametricTools::J<3>(smesh, eid).array() * GAUSSWEIGHTS<3>.array())).matrix() * initial_in_gp);
        }
    }

    // ******************** CG -> DG  ******************** //

  template <int CG, int DG>
  void CG2DG(const SasipMesh& smesh, CellVector<DG>& dg, const CGVector<CG>& cg)
    {
        // WHAT GAUSS DEGREE TO TAKE?
#define NGP 3

        assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == cg.rows());
        assert(static_cast<long int>(smesh.nx * smesh.ny) == dg.rows());

        const int cgshift = CG * smesh.nx + 1; //!< Index shift for each row

        // parallelize over the rows
#pragma omp parallel for
        for (size_t iy = 0; iy < smesh.ny; ++iy) {
            size_t dgi = smesh.nx * iy; //!< Index of dg vector
            size_t cgi = CG * cgshift * iy; //!< Lower left index of cg vector

            for (size_t ix = 0; ix < smesh.nx; ++ix, ++dgi, cgi += CG) {

	      Eigen::Matrix<double, (CG==2?9:4), 1> cg_local; //!< the 9 local unknowns in the element
	      if (CG==1)
		{
		  cg_local << cg(cgi), cg(cgi + 1),  cg(cgi + cgshift), cg(cgi + 1 + cgshift);
		}
	      else
		{
		  cg_local << cg(cgi), cg(cgi + 1), cg(cgi + 2), cg(cgi + cgshift), cg(cgi + 1 + cgshift),
                    cg(cgi + 2 + cgshift), cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift),
                    cg(cgi + 2 + 2 * cgshift);
		}
		
	      dg.row(dgi) = ParametricTools::massMatrix<DG>(smesh, dgi).inverse() * PSI<DG,NGP> * (ParametricTools::J<NGP>(smesh, dgi).array() * GAUSSWEIGHTS<NGP>.array() * (PHI<CG,NGP>.transpose() * cg_local).transpose().array()).matrix().transpose();
            }
        }
#undef NGP
    }

    // ******************** DG -> CG  ******************** //

  template<int DG>
    void DG2CGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<1>& cg_A, const CellVector<DG>& A)
    {
        const size_t CGDofsPerRow = smesh.nx + 1;
        const size_t cgi          = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)

	const Eigen::Matrix<double, 1, 4> At = A.row(c) * PSILagrange<DG,2>;
	
	cg_A(cgi) += 0.25 * At(0);
	cg_A(cgi + 1) += 0.25 * At(1);
	cg_A(cgi + CGDofsPerRow) += 0.25 * At(2);
	cg_A(cgi + CGDofsPerRow + 1) += 0.25 * At(3);
	
    }
    // void DG2CGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
    //     CGVector<1>& cg_A, const CellVector<3>& A)
    // {
    //     const size_t CGDofsPerRow = smesh.nx + 1;
    //     const size_t cgi
    //         = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)
    //     cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + 1) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + CGDofsPerRow) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    //     cg_A(cgi + CGDofsPerRow + 1) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    // }
    // void DG2CGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
    //     CGVector<1>& cg_A, const CellVector<6>& A)
    // {
    //     const size_t CGDofsPerRow = smesh.nx + 1;
    //     const size_t cgi          = CGDofsPerRow * cy + cx; //!< lower left index of CG-vector in element c = (cx,cy)

	
    //     cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + 1) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + CGDofsPerRow) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    //     cg_A(cgi + CGDofsPerRow + 1) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    // }
  template<int DG>
  void DG2CGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
        CGVector<2>& cg_A, const CellVector<DG>& A)
    {
        const size_t CGDofsPerRow = 2 * smesh.nx + 1;
        const size_t cgi = 2 * CGDofsPerRow * cy
            + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)

	const Eigen::Matrix<double, 1, 9> At = A.row(c) * PSILagrange<DG,3>;

	
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
    // void DG2CGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
    //     CGVector<2>& cg_A, const CellVector<3>& A)
    // {
    //     const size_t CGDofsPerRow = 2 * smesh.nx + 1;
    //     const size_t cgi = 2 * CGDofsPerRow * cy
    //         + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
    //     cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + 1) += 0.5 * (A(c, 0) - 0.5 * A(c, 2));
    //     cg_A(cgi + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + CGDofsPerRow) += 0.5 * (A(c, 0) - 0.5 * A(c, 1));
    //     cg_A(cgi + CGDofsPerRow + 1) += A(c, 0);
    //     cg_A(cgi + CGDofsPerRow + 2) += 0.5 * (A(c, 0) + 0.5 * A(c, 1));
    //     cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    //     cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * (A(c, 0) + 0.5 * A(c, 2));
    //     cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    // }
    // void DG2CGCell(const SasipMesh& smesh, const size_t c, const size_t cx, const size_t cy,
    //     CGVector<2>& cg_A, const CellVector<6>& A)
    // {
    //     const size_t CGDofsPerRow = 2 * smesh.nx + 1;
    //     const size_t cgi = 2 * CGDofsPerRow * cy
    //         + 2 * cx; //!< lower left index of CG-vector in element c = (cx,cy)
    //     cg_A(cgi) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + 1) += 0.5 * (A(c, 0) - 0.5 * A(c, 2));
    //     cg_A(cgi + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) - 0.5 * A(c, 2));
    //     cg_A(cgi + CGDofsPerRow) += 0.5 * (A(c, 0) - 0.5 * A(c, 1));
    //     cg_A(cgi + CGDofsPerRow + 1) += A(c, 0);
    //     cg_A(cgi + CGDofsPerRow + 2) += 0.5 * (A(c, 0) + 0.5 * A(c, 1));
    //     cg_A(cgi + 2 * CGDofsPerRow) += 0.25 * (A(c, 0) - 0.5 * A(c, 1) + 0.5 * A(c, 2));
    //     cg_A(cgi + 2 * CGDofsPerRow + 1) += 0.5 * (A(c, 0) + 0.5 * A(c, 2));
    //     cg_A(cgi + 2 * CGDofsPerRow + 2) += 0.25 * (A(c, 0) + 0.5 * A(c, 1) + 0.5 * A(c, 2));
    // }

    void DG2CGBoundary(const SasipMesh& smesh, CGVector<1>& cg_A)
    {
        const size_t CGDofsPerRow = smesh.nx + 1;
        const size_t UpperLeftIndex = CGDofsPerRow * smesh.ny;

#pragma omp parallel for
        for (size_t i = 0; i < smesh.nx + 1; ++i) {
            cg_A(i) *= 2.0;
            cg_A(UpperLeftIndex + i) *= 2.0;
        }
#pragma omp parallel for
        for (size_t i = 0; i < smesh.ny + 1; ++i) {
            cg_A(i * CGDofsPerRow) *= 2.0;
            cg_A(i * CGDofsPerRow + smesh.nx) *= 2.0;
        }
    }

    void DG2CGBoundary(const SasipMesh& smesh, CGVector<2>& cg_A)
    {
        const size_t CGDofsPerRow = 2 * smesh.nx + 1;
        const size_t UpperLeftIndex = 2 * CGDofsPerRow * smesh.ny;

#pragma omp parallel for
        for (size_t i = 0; i < 2 * smesh.nx + 1; ++i) {
            cg_A(i) *= 2.0;
            cg_A(UpperLeftIndex + i) *= 2.0;
        }
#pragma omp parallel for
        for (size_t i = 0; i < 2 * smesh.ny + 1; ++i) {
            cg_A(i * CGDofsPerRow) *= 2.0;
            cg_A(i * CGDofsPerRow + 2 * smesh.nx) *= 2.0;
        }
    }

    template <int CG, int DG>
    void DG2CG(const SasipMesh& smesh, CGVector<CG>& dest, const CellVector<DG>& src)
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

  template void DG2CG(const SasipMesh& smesh, CGVector<2>& dest, const CellVector<1>& src);
  template void DG2CG(const SasipMesh& smesh, CGVector<2>& dest, const CellVector<3>& src);
  template void DG2CG(const SasipMesh& smesh, CGVector<2>& dest, const CellVector<6>& src);
  template void DG2CG(const SasipMesh& smesh, CGVector<1>& dest, const CellVector<1>& src);
  template void DG2CG(const SasipMesh& smesh, CGVector<1>& dest, const CellVector<3>& src);
  template void DG2CG(const SasipMesh& smesh, CGVector<1>& dest, const CellVector<6>& src);
  
  template void CG2DG(const SasipMesh& smesh, CellVector<1>& dg, const CGVector<1>& cg);
  template void CG2DG(const SasipMesh& smesh, CellVector<3>& dg, const CGVector<1>& cg);
  template void CG2DG(const SasipMesh& smesh, CellVector<6>& dg, const CGVector<1>& cg);
  template void CG2DG(const SasipMesh& smesh, CellVector<1>& dg, const CGVector<2>& cg);
  template void CG2DG(const SasipMesh& smesh, CellVector<3>& dg, const CGVector<2>& cg);
  template void CG2DG(const SasipMesh& smesh, CellVector<6>& dg, const CGVector<2>& cg);
  

}
}
