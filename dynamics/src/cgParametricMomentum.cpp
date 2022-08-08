/*!
 * @file Mesh.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "cgParametricMomentum.hpp"
#include "ParametricTools.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "mevp.hpp"

#include "stopwatch.hpp"

namespace Nextsim {

extern Timer GlobalTimer;

// OLD MESH: Projection of a CG function to DG

template <>
template <>
void CGParametricMomentum<2, 8>::ProjectCGToDG(CellVector<3>& dg, const CGVector<2>& cg)
{
    // WHAT GAUSS DEGREE TO TAKE?

    assert(static_cast<long int>((2 * smesh.nx + 1) * (2 * smesh.ny + 1)) == cg.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == dg.rows());

    const int cgshift = 2 * smesh.nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t iy = 0; iy < smesh.ny; ++iy) {
        size_t dgi = smesh.nx * iy; //!< Index of dg vector
        size_t cgi = 2 * cgshift * iy; //!< Lower left index of cg vector

        for (size_t ix = 0; ix < smesh.nx; ++ix, ++dgi, cgi += 2) {

            Eigen::Matrix<double, 9, 1> cg_local; //!< the 9 local unknowns in the element
            cg_local << cg(cgi), cg(cgi + 1), cg(cgi + 2), cg(cgi + cgshift), cg(cgi + 1 + cgshift),
                cg(cgi + 2 + cgshift), cg(cgi + 2 * cgshift), cg(cgi + 1 + 2 * cgshift),
                cg(cgi + 2 + 2 * cgshift);

            dg.row(dgi) = ParametricTools::massMatrix<3>(smesh, dgi).inverse() * BiG33 * (ParametricTools::J<3>(smesh, dgi).array() * GAUSSWEIGHTS_3.array() * (CG_CG2FUNC_in_GAUSS3 * cg_local).transpose().array()).matrix().transpose();
        }
    }
}

////////////////////////////////////////////////// Strain (SasipMesh)

template <>
template <>
void CGParametricMomentum<2, 8>::ProjectCG2VelocityToDG1Strain(CellVector<8>& E11, CellVector<8>& E12, CellVector<8>& E22)
{
    assert(static_cast<long int>((2 * smesh.nx + 1) * (2 * smesh.ny + 1)) == vx.rows());
    assert(static_cast<long int>((2 * smesh.nx + 1) * (2 * smesh.ny + 1)) == vy.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E11.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E12.rows());
    assert(static_cast<long int>(smesh.nx * smesh.ny) == E22.rows());

    const int cgshift = 2 * smesh.nx + 1; //!< Index shift for each row

    if (precompute_matrices == 0) {
        // parallelize over the rows
#pragma omp parallel for
        for (size_t row = 0; row < smesh.ny; ++row) {
            int dgi = smesh.nx * row; //!< Index of dg vector
            int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

            for (size_t col = 0; col < smesh.nx; ++col, ++dgi, cgi += 2) {

                // get the 9 local x/y - velocity coefficients on the element
                Eigen::Matrix<double, 9, 1> vx_local;
                vx_local << vx(cgi), vx(cgi + 1), vx(cgi + 2), vx(cgi + cgshift), vx(cgi + 1 + cgshift),
                    vx(cgi + 2 + cgshift), vx(cgi + 2 * cgshift), vx(cgi + 1 + 2 * cgshift),
                    vx(cgi + 2 + 2 * cgshift);
                Eigen::Matrix<double, 9, 1> vy_local;
                vy_local << vy(cgi), vy(cgi + 1), vy(cgi + 2), vy(cgi + cgshift), vy(cgi + 1 + cgshift),
                    vy(cgi + 2 + cgshift), vy(cgi + 2 * cgshift), vy(cgi + 1 + 2 * cgshift),
                    vy(cgi + 2 + 2 * cgshift);

                // get the reference-Gradient of the velocity in the GP
                // (vx_i * d_x/y phi_i(q))
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> DX_VX_g = CG_CG2FUNC_DX_in_GAUSS3 * vx_local;
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> DY_VX_g = CG_CG2FUNC_DY_in_GAUSS3 * vx_local;
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> DX_VY_g = CG_CG2FUNC_DX_in_GAUSS3 * vy_local;
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> DY_VY_g = CG_CG2FUNC_DY_in_GAUSS3 * vy_local;

                const Eigen::Matrix<Nextsim::FloatType, 2, 9> dxT = (ParametricTools::dxT<3>(smesh, dgi).array().rowwise() * GAUSSWEIGHTS_3.array()).matrix();
                const Eigen::Matrix<Nextsim::FloatType, 2, 9> dyT = (ParametricTools::dyT<3>(smesh, dgi).array().rowwise() * GAUSSWEIGHTS_3.array()).matrix();

                // gradient of transformation
                //      [ dxT1, dyT1 ]     //            [ dyT2, -dxT2 ]
                // dT = 		       // J dT^{-T}=
                //      [ dxT2, dyT2 ]     //            [ -dyT1, dxT1 ]
                //

                const Eigen::Matrix<Nextsim::FloatType, 1, 9> J_dx_vx_g = dyT.row(1).array() * DX_VX_g.array() - dxT.row(1).array() * DY_VX_g.array();
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> J_dy_vx_g = dxT.row(0).array() * DY_VX_g.array() - dyT.row(0).array() * DX_VX_g.array();
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> J_dx_vy_g = dyT.row(1).array() * DX_VY_g.array() - dxT.row(1).array() * DY_VY_g.array();
                const Eigen::Matrix<Nextsim::FloatType, 1, 9> J_dy_vy_g = dxT.row(0).array() * DY_VY_g.array() - dyT.row(0).array() * DX_VY_g.array();

                const Eigen::Matrix<Nextsim::FloatType, 8, 8> imass = ParametricTools::massMatrix<8>(smesh, dgi).inverse();

                E11.row(dgi) = imass * (BiG83 * J_dx_vx_g.transpose());
                E22.row(dgi) = imass * (BiG83 * J_dy_vy_g.transpose());
                E12.row(dgi) = 0.5 * imass * (BiG83 * (J_dx_vy_g.transpose() + J_dy_vx_g.transpose()));
            }
        }
    } else if (precompute_matrices == 1) {
        // parallelize over the rows
#pragma omp parallel for
        for (size_t row = 0; row < smesh.ny; ++row) {
            int dgi = smesh.nx * row; //!< Index of dg vector
            int cgi = 2 * cgshift * row; //!< Lower left index of cg vector

            for (size_t col = 0; col < smesh.nx; ++col, ++dgi, cgi += 2) {

                // get the 9 local x/y - velocity coefficients on the element
                const Eigen::Matrix<double, 9, 1> vx_local({ vx(cgi), vx(cgi + 1), vx(cgi + 2), vx(cgi + cgshift), vx(cgi + 1 + cgshift),
                    vx(cgi + 2 + cgshift), vx(cgi + 2 * cgshift), vx(cgi + 1 + 2 * cgshift),
                    vx(cgi + 2 + 2 * cgshift) });
                const Eigen::Matrix<double, 9, 1> vy_local({ vy(cgi), vy(cgi + 1), vy(cgi + 2), vy(cgi + cgshift), vy(cgi + 1 + cgshift),
                    vy(cgi + 2 + cgshift), vy(cgi + 2 * cgshift), vy(cgi + 1 + 2 * cgshift),
                    vy(cgi + 2 + 2 * cgshift) });

                E11.row(dgi) = ptrans.iMgradX[dgi] * vx_local;
                E22.row(dgi) = ptrans.iMgradY[dgi] * vy_local;
                E12.row(dgi) = 0.5 * (ptrans.iMgradX[dgi] * vy_local + ptrans.iMgradY[dgi] * vx_local);
            }
        }
    }
}

////////////////////////////////////////////////// STRESS Tensor
// Sasip-Mesh Interface
template <int CG, int DGstress>
template <int DG>
void CGParametricMomentum<CG, DGstress>::AddStressTensor(const double scale, CGVector<CG>& tx,
    CGVector<CG>& ty, const CellVector<DG>& S11, const CellVector<DG>& S12,
    const CellVector<DG>& S22) const
{
    // parallelization in tripes
    for (size_t p = 0; p < 2; ++p)
#pragma omp parallel for schedule(static)
        for (size_t cy = 0; cy < smesh.ny; ++cy) //!< loop over all cells of the mesh
        {
            if (cy % 2 == p) {
                size_t c = smesh.nx * cy;
                for (size_t cx = 0; cx < smesh.nx; ++cx, ++c) //!< loop over all cells of the mesh
                    AddStressTensorCell(scale, c, cx, cy, tx, ty, S11, S12, S22);
            }
        }
}

//! Sets the vector to zero along the boundary
template <>
void CGParametricMomentum<1, 3>::DirichletZero(CGVector<1>& v)
{

    size_t upperleftindex = (smesh.nx + 1) * smesh.ny;
    for (size_t i = 0; i < smesh.nx + 1; ++i) {
        v(i, 0) = 0.0;
        v(upperleftindex + i, 0) = 0.0;
    }
    size_t indecesperrow = smesh.nx + 1;
    size_t lowerrightindex = smesh.nx;
    for (size_t i = 0; i < smesh.ny + 1; ++i) {
        v(indecesperrow * i, 0) = 0.0;
        v(lowerrightindex + indecesperrow * i, 0) = 0.0;
    }
}
template <>
void CGParametricMomentum<2, 8>::DirichletZero(CGVector<2>& v)
{

    size_t upperleftindex = (2 * smesh.nx + 1) * 2 * smesh.ny;
    for (size_t i = 0; i < 2 * smesh.nx + 1; ++i) {
        v(i, 0) = 0.0;
        v(upperleftindex + i, 0) = 0.0;
    }
    size_t indecesperrow = 2 * smesh.nx + 1;
    size_t lowerleftindex = 2 * smesh.nx;
    for (size_t i = 0; i < 2 * smesh.ny + 1; ++i) {
        v(indecesperrow * i, 0) = 0.0;
        v(lowerleftindex + indecesperrow * i, 0) = 0.0;
    }
}

template <int CG, int DGstress>
template <int DG>
void CGParametricMomentum<CG, DGstress>::InterpolateDGToCG(
    CGVector<CG>& cg_A, const CellVector<DG>& A) const
{
    cg_A.zero();

    // parallelization by running over stripes
    for (size_t p = 0; p < 2; ++p) {
#pragma omp parallel for
        for (size_t cy = 0; cy < smesh.ny; ++cy) {
            if (cy % 2 == p)
                continue;

            size_t c = cy * smesh.nx;

            for (size_t cx = 0; cx < smesh.nx; ++cx, ++c)
                InterpolateDGToCGCell(c, cx, cy, cg_A, A);
        }
    }

    // InterpolateDGToCGCell adds to the CG-Dofs with correct weighting. Then we adjust the boundary
    InterpolateDGToCGBoundary(cg_A);
}

// --------------------------------------------------

template <int CG, int DGstress>
template <int DG>
void CGParametricMomentum<CG, DGstress>::mEVPIteration(const size_t NT_evp, const double alpha, const double beta,
    double dt_adv,
    const CellVector<DG>& H, const CellVector<DG>& A,
    CellVector<DGstress>& E11, CellVector<DGstress>& E12, CellVector<DGstress>& E22,
    CellVector<DGstress>& S11, CellVector<DGstress>& S12, CellVector<DGstress>& S22)
{

    // copy old velocity
    CGVector<CG> vx_mevp = vx;
    CGVector<CG> vy_mevp = vy;

    // interpolate ice height and concentration to local cg Variables
    InterpolateDGToCG(cg_A, A); // should not be within this class
    InterpolateDGToCG(cg_H, H);
    // limit A to [0,1] and H to [0, ...)
    cg_A = cg_A.cwiseMin(1.0);
    cg_A = cg_A.cwiseMax(0.0);
    cg_H = cg_H.cwiseMax(1.e-4);

    // MEVP subcycling
    for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {

        Nextsim::GlobalTimer.start("time loop - mevp - strain");
        //! Compute Strain Rate
        // momentum.ProjectCG2VelocityToDG1Strain(ptrans_stress, E11, E12, E22);
        ProjectCG2VelocityToDG1Strain(E11, E12, E22);
        Nextsim::GlobalTimer.stop("time loop - mevp - strain");

        Nextsim::GlobalTimer.start("time loop - mevp - stress");
        Nextsim::mEVP::StressUpdateHighOrder(smesh, S11, S12, S22, E11, E12, E22, H, A, alpha, beta);
        // Nextsim::mEVP::StressUpdateHighOrder(ptrans_stress, smesh, S11, S12, S22, E11, E12, E22, H, A,
        //     EVPParameters::Pstar, EVPParameters::DeltaMin, alpha, beta);
        Nextsim::GlobalTimer.stop("time loop - mevp - stress");

        Nextsim::GlobalTimer.start("time loop - mevp - update");
        //! Update
        Nextsim::GlobalTimer.start("time loop - mevp - update1");

        //	    update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
        for (int i = 0; i < vx.rows(); ++i) {
            vx(i) = (1.0
                / (EVPParameters::rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A(i) * EVPParameters::F_ocean
                        * fabs(ox(i) - vx(i))) // implicit parts
                * (EVPParameters::rho_ice * cg_H(i) / dt_adv
                        * (beta * vx(i) + vx_mevp(i))
                    + // pseudo-timestepping
                    cg_A(i)
                        * (EVPParameters::F_atm * fabs(ax(i)) * ax(i) + // atm forcing
                            EVPParameters::F_ocean * fabs(ox(i) - vx(i))
                                * ox(i)) // ocean forcing
                    + EVPParameters::rho_ice * cg_H(i) * EVPParameters::fc
                        * (vy(i) - oy(i)) // cor + surface
                    ));
            vy(i) = (1.0
                / (EVPParameters::rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A(i) * EVPParameters::F_ocean
                        * fabs(oy(i) - vy(i))) // implicit parts
                * (EVPParameters::rho_ice * cg_H(i) / dt_adv
                        * (beta * vy(i) + vy_mevp(i))
                    + // pseudo-timestepping
                    cg_A(i)
                        * (EVPParameters::F_atm * fabs(ay(i)) * ay(i) + // atm forcing
                            EVPParameters::F_ocean * fabs(oy(i) - vy(i))
                                * oy(i)) // ocean forcing
                    + EVPParameters::rho_ice * cg_H(i) * EVPParameters::fc
                        * (ox(i) - vx(i))));
        }
        Nextsim::GlobalTimer.stop("time loop - mevp - update1");

        Nextsim::GlobalTimer.start("time loop - mevp - update2");
        // Implicit etwas ineffizient
#pragma omp parallel for
        for (int i = 0; i < tmpx.rows(); ++i)
            tmpx(i) = tmpy(i) = 0;

        Nextsim::GlobalTimer.start("time loop - mevp - update2 -stress");
        // AddStressTensor(ptrans_stress, -1.0, tmpx, tmpy, S11, S12, S22);
        AddStressTensor(-1.0, tmpx, tmpy, S11, S12, S22);
        Nextsim::GlobalTimer.stop("time loop - mevp - update2 -stress");

#pragma omp parallel for
        for (int i = 0; i < vx.rows(); ++i) {
            vx(i) += (1.0
                         / (EVPParameters::rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                             + cg_A(i) * EVPParameters::F_ocean
                                 * fabs(ox(i) - vx(i))) // implicit parts
                         * tmpx(i))
                / lumpedcgmass(i);
            ;

            vy(i) += (1.0
                         / (EVPParameters::rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                             + cg_A(i) * EVPParameters::F_ocean
                                 * fabs(oy(i) - vy(i))) // implicit parts
                         * tmpy(i))
                / lumpedcgmass(i);
            ;
        }
        Nextsim::GlobalTimer.stop("time loop - mevp - update2");
        Nextsim::GlobalTimer.stop("time loop - mevp - update");

        Nextsim::GlobalTimer.start("time loop - mevp - bound.");
        DirichletZero();
        Nextsim::GlobalTimer.stop("time loop - mevp - bound.");
    }
}

// --------------------------------------------------

// template class CGParametricMomentum<1, 3>;
template class CGParametricMomentum<2, 8>;

template void CGParametricMomentum<2, 8>::AddStressTensor(const double scale, CGVector<2>& tx,
    CGVector<2>& ty, const CellVector<8>& S11, const CellVector<8>& S12,
    const CellVector<8>& S22) const;

// --------------------------------------------------

template void CGParametricMomentum<2, 8>::mEVPIteration(size_t NT_evp, double alpha, double beta,
    double dt_adv,
    const CellVector<3>& H, const CellVector<3>& A,
    CellVector<8>& E11, CellVector<8>& E12, CellVector<8>& E22,
    CellVector<8>& S11, CellVector<8>& S12, CellVector<8>& S22);

// template void CGParametricMomentum<1, 3>::InterpolateDGToCG(
//     CGVector<1>& cg_A, const CellVector<1>& A) const;
// template void CGParametricMomentum<1, 3>::InterpolateDGToCG(
//     CGVector<1>& cg_A, const CellVector<3>& A) const;
template void CGParametricMomentum<2, 8>::InterpolateDGToCG(
    CGVector<2>& cg_A, const CellVector<1>& A) const;
template void CGParametricMomentum<2, 8>::InterpolateDGToCG(
    CGVector<2>& cg_A, const CellVector<3>& A) const;
template void CGParametricMomentum<2, 8>::InterpolateDGToCG(
    CGVector<2>& cg_A, const CellVector<6>& A) const;

} /* namespace Nextsim */
