/*!
 * @file ParametricMomentum.cpp
 * @date 06 Dec 2024
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "cgParametricMomentum.hpp"
#include "Interpolations.hpp"
#include "ParametricTools.hpp"
#include "VectorManipulations.hpp"
#include "dgVisu.hpp"

namespace Nextsim {

#define DGSTRESS(CG) ((CG == 1 ? 3 : (CG == 2 ? 8 : -1)))

////////////////////////////////////////////////// Strain (ParametricMesh)

template <int CG> void CGParametricMomentum<CG>::ProjectCGVelocityToDGStrain()
{
    // !!! must still be converted to the spherical system!!!

    assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == vx.rows());
    assert(static_cast<long int>((CG * smesh.nx + 1) * (CG * smesh.ny + 1)) == vy.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E11.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E12.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E22.rows());

    const int cgshift = CG * smesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < smesh.ny; ++row) {
        int dgi = smesh.nx * row; //!< Index of dg vector
        int cgi = CG * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < smesh.nx; ++col, ++dgi, cgi += CG) { // loop over all elements

            if (smesh.landmask[dgi] == 0) // only on ice
                continue;

            // get the 4 (cg1) 9 (cg2) local x/y - velocity coefficients on the element
            Eigen::Matrix<double, CGDOFS(CG), 1> vx_local, vy_local;
            if (CG == 1) {
                vx_local << vx(cgi), vx(cgi + 1), vx(cgi + cgshift), vx(cgi + 1 + cgshift);
                vy_local << vy(cgi), vy(cgi + 1), vy(cgi + cgshift), vy(cgi + 1 + cgshift);
            } else if (CG == 2) {
                vx_local << vx(cgi), vx(cgi + 1), vx(cgi + 2), vx(cgi + cgshift),
                    vx(cgi + 1 + cgshift), vx(cgi + 2 + cgshift), vx(cgi + 2 * cgshift),
                    vx(cgi + 1 + 2 * cgshift), vx(cgi + 2 + 2 * cgshift);

                vy_local << vy(cgi), vy(cgi + 1), vy(cgi + 2), vy(cgi + cgshift),
                    vy(cgi + 1 + cgshift), vy(cgi + 2 + cgshift), vy(cgi + 2 * cgshift),
                    vy(cgi + 1 + 2 * cgshift), vy(cgi + 2 + 2 * cgshift);
            } else
                abort();

            // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
            // by integrating rhs and inverting with dG(stress) mass matrix
            //
            E11.row(dgi) = pmap.iMgradX[dgi] * vx_local;
            E22.row(dgi) = pmap.iMgradY[dgi] * vy_local;
            E12.row(dgi) = 0.5 * (pmap.iMgradX[dgi] * vy_local + pmap.iMgradY[dgi] * vx_local);

            if (smesh.CoordinateSystem == SPHERICAL) {
                E11.row(dgi) -= pmap.iMM[dgi] * vy_local;
                E12.row(dgi) += 0.5 * pmap.iMM[dgi] * vx_local;
            }
        }
    }
}

////////////////////////////////////////////////// STRESS Tensor
// Sasip-Mesh Interface
template <int CG>
void CGParametricMomentum<CG>::DivergenceOfStress(
    const double scale, CGVector<CG>& tx, CGVector<CG>& ty) const
{
#pragma omp parallel for
    for (size_t i = 0; i < tx.rows(); ++i) {
        tx(i) = 0.0;
        ty(i) = 0.0;
    }

    // parallelization in stripes
    for (size_t p = 0; p < 2; ++p)
#pragma omp parallel for schedule(static)
        for (size_t cy = 0; cy < smesh.ny; ++cy) //!< loop over all cells of the mesh
        {
            if (cy % 2 == p) {
                size_t c = smesh.nx * cy;
                for (size_t cx = 0; cx < smesh.nx; ++cx, ++c) //!< loop over all cells of the mesh
                    if (smesh.landmask[c] == 1) // only on ice!
                        AddStressTensorCell(scale, c, cx, cy, tx, ty);
            }
        }
    // set zero on the Dirichlet boundaries
    DirichletZero(tx);
    DirichletZero(ty);
    // add the contributions on the periodic boundaries
    VectorManipulations::CGAveragePeriodic(smesh, tx);
    VectorManipulations::CGAveragePeriodic(smesh, ty);
}

template <int CG> void CGParametricMomentum<CG>::DirichletZero(CGVector<CG>& v) const
{
    // the four segments bottom, right, top, left, are each processed in parallel
    for (size_t seg = 0; seg < 4; ++seg) {
#pragma omp parallel for
        for (size_t i = 0; i < smesh.dirichlet[seg].size(); ++i) {

            const size_t eid = smesh.dirichlet[seg][i];
            const size_t ix = eid % smesh.nx; // compute 'coordinate' of element
            const size_t iy = eid / smesh.nx;

            if (seg == 0) // bottom
                for (size_t j = 0; j < CG + 1; ++j)
                    v(iy * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) = 0.0;
            else if (seg == 1) // right
                for (size_t j = 0; j < CG + 1; ++j)
                    v(iy * CG * (CG * smesh.nx + 1) + CG * ix + CG + (CG * smesh.nx + 1) * j, 0)
                        = 0.0;
            else if (seg == 2) // top
                for (size_t j = 0; j < CG + 1; ++j)
                    v((iy + 1) * CG * (CG * smesh.nx + 1) + CG * ix + j, 0) = 0.0;
            else if (seg == 3) // left
                for (size_t j = 0; j < CG + 1; ++j)
                    v(iy * CG * (CG * smesh.nx + 1) + CG * ix + (CG * smesh.nx + 1) * j, 0) = 0.0;
            else {
                std::cerr << "That should not have happened!" << std::endl;
                abort();
            }
        }
    }
}

template <int CG> void CGParametricMomentum<CG>::CheckPeriodicity(CGVector<CG>& v)
{
    // the two segments bottom, right, top, left, are each processed in parallel
    for (size_t seg = 0; seg < smesh.periodic.size(); ++seg) {
        // #pragma omp parallel for
        for (size_t i = 0; i < smesh.periodic[seg].size(); ++i) {

            const size_t ptype = smesh.periodic[seg][i][0];
            const size_t eid_lb = smesh.periodic[seg][i][2];
            const size_t eid_rt = smesh.periodic[seg][i][1];

            size_t ix_lb = eid_lb % smesh.nx;
            size_t iy_lb = eid_lb / smesh.nx;
            size_t i0_lb = (CG * smesh.nx + 1) * CG * iy_lb
                + CG * ix_lb; // lower/left index in left/bottom element
            size_t ix_rt = eid_rt % smesh.nx;
            size_t iy_rt = eid_rt / smesh.nx;
            size_t i0_rt = (CG * smesh.nx + 1) * CG * iy_rt
                + CG * ix_rt; // lower/left index in right/top element

            std::cout << std::setprecision(16);
            if (ptype == 0) // X-edge, bottom/top
            {
                for (size_t j = 0; j <= CG; ++j) {
                    double check = v(i0_lb + j) - v(i0_rt + CG * (CG * smesh.nx + 1) + j);
                    if (fabs(check) > 1.e-13)
                        std::cout << v(i0_lb + j) - v(i0_rt + CG * (CG * smesh.nx + 1) + j)
                                  << std::endl;
                }
            } else if (ptype == 1) // Y-edge, left/right
            {
                for (size_t j = 0; j < CG; ++j) {
                    double check = v(i0_lb + j * (CG * smesh.nx + 1))
                        - v(i0_rt + CG + j * (CG * smesh.nx + 1));
                    if (fabs(check) > 1.e-13)
                        std::cout << v(i0_lb + j * (CG * smesh.nx + 1))
                                - v(i0_rt + CG + j * (CG * smesh.nx + 1))
                                  << std::endl;
                }
            } else
                abort();
        }
    }
}

// --------------------------------------------------

template <int CG>
template <int DG>
void CGParametricMomentum<CG>::prepareIteration(const DGVector<DG>& H, const DGVector<DG>& A)
{
    // copy old velocity
    vx_mevp = vx;
    vy_mevp = vy;
    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(smesh, cg_A, A);
    VectorManipulations::CGAveragePeriodic(smesh, cg_A);
    Interpolations::DG2CG(smesh, cg_H, H);
    VectorManipulations::CGAveragePeriodic(smesh, cg_H);

    // limit A to [0,1] and H to [0, ...)
    cg_A = cg_A.cwiseMin(1.0);
    cg_A = cg_A.cwiseMax(1.e-4);
    cg_H = cg_H.cwiseMax(1.e-4);
}

// --------------------------------------------------
template <int CG>
template <int DG>
void CGParametricMomentum<CG>::prepareIteration(
    const DGVector<DG>& H, const DGVector<DG>& A, const DGVector<DG>& D)
{

    // set the average sub-iteration velocity to zero
    avg_vx.setZero();
    avg_vy.setZero();

    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(smesh, cg_A, A);
    VectorManipulations::CGAveragePeriodic(smesh, cg_A);
    Interpolations::DG2CG(smesh, cg_H, H);
    VectorManipulations::CGAveragePeriodic(smesh, cg_H);
    Interpolations::DG2CG(smesh, cg_D, D);
    VectorManipulations::CGAveragePeriodic(smesh, cg_D);

    // limit A and D to [0,1] and H to [0, ...)
    cg_A = cg_A.cwiseMin(1.0);
    cg_A = cg_A.cwiseMax(1.e-4);
    cg_H = cg_H.cwiseMax(1.e-4);
    cg_D = cg_D.cwiseMin(1.0);
    cg_D = cg_D.cwiseMax(0.0);
}

// --------------------------------------------------

template class CGParametricMomentum<1>;
template class CGParametricMomentum<2>;

template void CGParametricMomentum<1>::prepareIteration(const DGVector<1>& H, const DGVector<1>& A);
template void CGParametricMomentum<1>::prepareIteration(const DGVector<3>& H, const DGVector<3>& A);
template void CGParametricMomentum<1>::prepareIteration(const DGVector<6>& H, const DGVector<6>& A);
template void CGParametricMomentum<2>::prepareIteration(const DGVector<1>& H, const DGVector<1>& A);
template void CGParametricMomentum<2>::prepareIteration(const DGVector<3>& H, const DGVector<3>& A);
template void CGParametricMomentum<2>::prepareIteration(const DGVector<6>& H, const DGVector<6>& A);

template void CGParametricMomentum<1>::prepareIteration(
    const DGVector<1>& H, const DGVector<1>& A, const DGVector<1>& D);
template void CGParametricMomentum<1>::prepareIteration(
    const DGVector<3>& H, const DGVector<3>& A, const DGVector<3>& D);
template void CGParametricMomentum<1>::prepareIteration(
    const DGVector<6>& H, const DGVector<6>& A, const DGVector<6>& D);
template void CGParametricMomentum<2>::prepareIteration(
    const DGVector<1>& H, const DGVector<1>& A, const DGVector<1>& D);
template void CGParametricMomentum<2>::prepareIteration(
    const DGVector<3>& H, const DGVector<3>& A, const DGVector<3>& D);
template void CGParametricMomentum<2>::prepareIteration(
    const DGVector<6>& H, const DGVector<6>& A, const DGVector<6>& D);

} /* namespace Nextsim */
