/*!
 * @file Mesh.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "cgParametricMomentum.hpp"
#include "ParametricTools.hpp"
#include "Interpolations.hpp"
#include "codeGenerationCGinGauss.hpp"
#include "mevp.hpp"

#include "stopwatch.hpp"

namespace Nextsim {

extern Timer GlobalTimer;


////////////////////////////////////////////////// Strain (SasipMesh)

template <>
void CGParametricMomentum<2, 8>::ProjectCG2VelocityToDG1Strain()
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
void CGParametricMomentum<CG, DGstress>::AddStressTensor(const double scale, CGVector<CG>& tx,
    CGVector<CG>& ty) const
{
    // parallelization in tripes
    for (size_t p = 0; p < 2; ++p)
#pragma omp parallel for schedule(static)
        for (size_t cy = 0; cy < smesh.ny; ++cy) //!< loop over all cells of the mesh
        {
            if (cy % 2 == p) {
                size_t c = smesh.nx * cy;
                for (size_t cx = 0; cx < smesh.nx; ++cx, ++c) //!< loop over all cells of the mesh
		  AddStressTensorCell(scale, c, cx, cy, tx, ty);
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

// --------------------------------------------------

template <int CG, int DGstress>
template <int DG>
void CGParametricMomentum<CG, DGstress>::mEVPIteration(const VPParameters& vpparameters,
						       const size_t NT_evp, const double alpha, const double beta,
    double dt_adv,
						       const CellVector<DG>& H, const CellVector<DG>& A)
{

    // copy old velocity
    CGVector<CG> vx_mevp = vx;
    CGVector<CG> vy_mevp = vy;

    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(smesh, cg_A, A);
    Interpolations::DG2CG(smesh, cg_H, H);
    // limit A to [0,1] and H to [0, ...)
    cg_A = cg_A.cwiseMin(1.0);
    cg_A = cg_A.cwiseMax(0.0);
    cg_H = cg_H.cwiseMax(1.e-4);

    // MEVP subcycling
    for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {

        Nextsim::GlobalTimer.start("time loop - mevp - strain");
        //! Compute Strain Rate
        // momentum.ProjectCG2VelocityToDG1Strain(ptrans_stress, E11, E12, E22);
        ProjectCG2VelocityToDG1Strain();
        Nextsim::GlobalTimer.stop("time loop - mevp - strain");

        Nextsim::GlobalTimer.start("time loop - mevp - stress");
        Nextsim::mEVP::StressUpdateHighOrder(vpparameters,smesh, S11, S12, S22, E11, E12, E22, H, A, alpha, beta);
        // Nextsim::mEVP::StressUpdateHighOrder(ptrans_stress, smesh, S11, S12, S22, E11, E12, E22, H, A,
        //     vpparameters.Pstar, vpparameters.DeltaMin, alpha, beta);
        Nextsim::GlobalTimer.stop("time loop - mevp - stress");

        Nextsim::GlobalTimer.start("time loop - mevp - update");
        //! Update
        Nextsim::GlobalTimer.start("time loop - mevp - update1");

        //	    update by a loop.. implicit parts and h-dependent
#pragma omp parallel for
        for (int i = 0; i < vx.rows(); ++i) {
            vx(i) = (1.0
                / (vpparameters.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A(i) * vpparameters.F_ocean
                        * fabs(ox(i) - vx(i))) // implicit parts
                * (vpparameters.rho_ice * cg_H(i) / dt_adv
                        * (beta * vx(i) + vx_mevp(i))
                    + // pseudo-timestepping
                    cg_A(i)
                        * (vpparameters.F_atm * fabs(ax(i)) * ax(i) + // atm forcing
                            vpparameters.F_ocean * fabs(ox(i) - vx(i))
                                * ox(i)) // ocean forcing
                    + vpparameters.rho_ice * cg_H(i) * vpparameters.fc
                        * (vy(i) - oy(i)) // cor + surface
                    ));
            vy(i) = (1.0
                / (vpparameters.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                    + cg_A(i) * vpparameters.F_ocean
                        * fabs(oy(i) - vy(i))) // implicit parts
                * (vpparameters.rho_ice * cg_H(i) / dt_adv
                        * (beta * vy(i) + vy_mevp(i))
                    + // pseudo-timestepping
                    cg_A(i)
                        * (vpparameters.F_atm * fabs(ay(i)) * ay(i) + // atm forcing
                            vpparameters.F_ocean * fabs(oy(i) - vy(i))
                                * oy(i)) // ocean forcing
                    + vpparameters.rho_ice * cg_H(i) * vpparameters.fc
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
        AddStressTensor(-1.0, tmpx, tmpy);
        Nextsim::GlobalTimer.stop("time loop - mevp - update2 -stress");

#pragma omp parallel for
        for (int i = 0; i < vx.rows(); ++i) {
            vx(i) += (1.0
                         / (vpparameters.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                             + cg_A(i) * vpparameters.F_ocean
                                 * fabs(ox(i) - vx(i))) // implicit parts
                         * tmpx(i))
                / lumpedcgmass(i);
            ;

            vy(i) += (1.0
                         / (vpparameters.rho_ice * cg_H(i) / dt_adv * (1.0 + beta) // implicit parts
                             + cg_A(i) * vpparameters.F_ocean
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


// --------------------------------------------------

template void CGParametricMomentum<2, 8>::mEVPIteration(const VPParameters& vpparameters,
							size_t NT_evp, double alpha, double beta,
    double dt_adv,
							const CellVector<3>& H, const CellVector<3>& A);


} /* namespace Nextsim */
