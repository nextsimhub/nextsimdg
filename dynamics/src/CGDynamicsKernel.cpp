/*!
 * @file CGDynamicsKernel.cpp
 *
 * @date 27 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar ÓLason <einar.olason@nersc.no>
 */

/*
 * The implementation of DynamicsKernel which uses continuous Galerkin (CG) numerics.
 */

#include "include/CGDynamicsKernel.hpp"
#include "include/ModelArray.hpp"

#include "include/Interpolations.hpp"
#include "include/ParametricMap.hpp"
#include "include/VectorManipulations.hpp"
#include "include/cgVector.hpp"

namespace Nextsim {

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    DynamicsKernel<DGadvection, DGstressComp>::initialise(coords, isSpherical, mask);

    //! Initialize the parametric momentum map
    pmap = new ParametricMomentumMap<CGdegree>(*smesh);
    pmap->InitializeLumpedCGMassMatrix();
    pmap->InitializeDivSMatrices();

    u.resize_by_mesh(*smesh);
    v.resize_by_mesh(*smesh);

    cgH.resize_by_mesh(*smesh);
    cgA.resize_by_mesh(*smesh);

    uGradSeasurfaceHeight.resize_by_mesh(*smesh);
    vGradSeasurfaceHeight.resize_by_mesh(*smesh);

    dStressX.resize_by_mesh(*smesh);
    dStressY.resize_by_mesh(*smesh);

    uOcean.resize_by_mesh(*smesh);
    vOcean.resize_by_mesh(*smesh);

    uAtmos.resize_by_mesh(*smesh);
    vAtmos.resize_by_mesh(*smesh);
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::setData(const std::string& name, const ModelArray& data)
{
    if (name == uName) {
        // FIXME take into account possibility to restart form CG
        // CGModelArray::ma2cg(data, u);
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, u, utmp);
    } else if (name == vName) {
        // CGModelArray::ma2cg(data, v);
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, v, vtmp);
    } else if (name == uWindName) {
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, uAtmos, utmp);
    } else if (name == vWindName) {
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, vAtmos, vtmp);
    } else if (name == uOceanName) {
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, uOcean, utmp);
    } else if (name == vOceanName) {
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, vOcean, vtmp);
    } else {
        DynamicsKernel<DGadvection, DGstressComp>::setData(name, data);
    }
}

template <int DGadvection>
ModelArray CGDynamicsKernel<DGadvection>::getDG0Data(const std::string& name) const
{
    if (name == uName) {
        ModelArray data(ModelArray::Type::U);
        DGVector<DGadvection> utmp(*smesh);
        Nextsim::Interpolations::CG2DG(*smesh, utmp, u);
        return DGModelArray::dg2ma(utmp, data);
    } else if (name == vName) {
        ModelArray data(ModelArray::Type::V);
        DGVector<DGadvection> vtmp(*smesh);
        Nextsim::Interpolations::CG2DG(*smesh, vtmp, v);
        return DGModelArray::dg2ma(vtmp, data);
    } else if (name == uIOStressName) {
        ModelArray data(ModelArray::Type::U);
        DGVector<DGadvection> utmp(*smesh);
        Nextsim::Interpolations::CG2DG(*smesh, utmp, getIceOceanStress().first);
        return DGModelArray::dg2ma(utmp, data);
    } else if (name == vIOStressName) {
        ModelArray data(ModelArray::Type::V);
        DGVector<DGadvection> vtmp(*smesh);
        Nextsim::Interpolations::CG2DG(*smesh, vtmp, getIceOceanStress().second);
        return DGModelArray::dg2ma(vtmp, data);
    } else {
        return DynamicsKernel<DGadvection, DGstressComp>::getDG0Data(name);
    }
}

template <int DGadvection> void CGDynamicsKernel<DGadvection>::prepareAdvection()
{
    dgtransport->prepareAdvection(u, v);
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::ComputeGradientOfSeaSurfaceHeight(
    const DGVector<1>& seaSurfaceHeight)
{
    // First transform to CG1 Vector and do all computations in CG1
    CGVector<1> cgSeasurfaceHeight(*smesh);
    Interpolations::DG2CG(*smesh, cgSeasurfaceHeight, seaSurfaceHeight);

    CGVector<1> uGrad(*smesh);
    CGVector<1> vGrad(*smesh);
    uGrad.setZero();
    vGrad.setZero();

    // parallelization in stripes
    for (size_t p = 0; p < 2; ++p) {
#pragma omp parallel for schedule(static)
        for (size_t cy = 0; cy < smesh->ny; ++cy) {
            //!< loop over all cells of the mesh
            if (cy % 2 == p) {
                size_t eid = smesh->nx * cy;
                size_t cg1id = cy * (smesh->nx + 1);
                for (size_t cx = 0; cx < smesh->nx; ++cx, ++eid, ++cg1id) {
                    // get local CG nodes
                    Eigen::Vector<Nextsim::FloatType, 4> loc_cgSSH = { cgSeasurfaceHeight(cg1id),
                        cgSeasurfaceHeight(cg1id + 1), cgSeasurfaceHeight(cg1id + smesh->nx + 1),
                        cgSeasurfaceHeight(cg1id + smesh->nx + 1 + 1) };

                    // compute grad
                    Eigen::Vector<Nextsim::FloatType, 4> tx = pmap->dX_SSH[eid] * loc_cgSSH;
                    Eigen::Vector<Nextsim::FloatType, 4> ty = pmap->dY_SSH[eid] * loc_cgSSH;

                    // add global vector
                    uGrad(cg1id) -= tx(0);
                    uGrad(cg1id + 1) -= tx(1);
                    uGrad(cg1id + smesh->nx + 1) -= tx(2);
                    uGrad(cg1id + smesh->nx + 1 + 1) -= tx(3);
                    vGrad(cg1id) -= ty(0);
                    vGrad(cg1id + 1) -= ty(1);
                    vGrad(cg1id + smesh->nx + 1) -= ty(2);
                    vGrad(cg1id + smesh->nx + 1 + 1) -= ty(3);
                }
            }
        }
    }

    // scale with mass
#pragma omp parallel for
    for (size_t i = 0; i < uGrad.rows(); ++i) {
        uGrad(i) /= pmap->lumpedcg1mass(i);
        vGrad(i) /= pmap->lumpedcg1mass(i);
    }

    ///// correct boundary (just extend in last elements)
    size_t cg1row = smesh->nx + 1;
    size_t topleft = smesh->ny * cg1row;
    for (size_t i = 1; i < smesh->nx; ++i) // bottom / top
    {
        uGrad(i) = uGrad(i + cg1row);
        vGrad(i) = vGrad(i + cg1row);
        uGrad(topleft + i) = uGrad(topleft + i - cg1row);
        vGrad(topleft + i) = vGrad(topleft + i - cg1row);
    }
    for (size_t i = 1; i < smesh->ny; ++i) // left / right
    {
        uGrad(i * cg1row) = uGrad(i * cg1row + 1);
        vGrad(i * cg1row) = vGrad(i * cg1row + 1);

        uGrad(i * cg1row + cg1row - 1) = uGrad(i * cg1row + cg1row - 1 - 1);
        vGrad(i * cg1row + cg1row - 1) = vGrad(i * cg1row + cg1row - 1 - 1);
    }
    // corners
    uGrad(0) = uGrad(cg1row + 1);
    vGrad(0) = vGrad(cg1row + 1);
    uGrad(smesh->nx) = uGrad(smesh->nx + cg1row - 1);
    vGrad(smesh->nx) = vGrad(smesh->nx + cg1row - 1);
    uGrad(smesh->ny * cg1row) = uGrad((smesh->ny - 1) * cg1row + 1);
    vGrad(smesh->ny * cg1row) = vGrad((smesh->ny - 1) * cg1row + 1);
    uGrad((smesh->ny + 1) * cg1row - 1) = uGrad((smesh->ny) * cg1row - 1 - 1);
    vGrad((smesh->ny + 1) * cg1row - 1) = vGrad((smesh->ny) * cg1row - 1 - 1);

    // Interpolate to CG2 (maybe own function in interpolation?)
    if (CGDEGREE == 1) {
        uGradSeasurfaceHeight = uGrad;
        vGradSeasurfaceHeight = vGrad;
    } else {
        // outer nodes
        size_t icg1 = 0;
#pragma omp parallel for
        for (size_t iy = 0; iy <= smesh->ny; ++iy) {
            size_t icg2 = (2 * smesh->nx + 1) * 2 * iy;
            for (size_t ix = 0; ix <= smesh->nx; ++ix, ++icg1, icg2 += 2) {
                uGradSeasurfaceHeight(icg2) = uGrad(icg1);
                vGradSeasurfaceHeight(icg2) = vGrad(icg1);
            }
        }
        // along lines
#pragma omp parallel for
        for (size_t iy = 0; iy <= smesh->ny; ++iy) // horizontal
        {
            size_t icg1 = (smesh->nx + 1) * iy;
            size_t icg2 = (2 * smesh->nx + 1) * 2 * iy + 1;
            for (size_t ix = 0; ix < smesh->nx; ++ix, ++icg1, icg2 += 2) {
                uGradSeasurfaceHeight(icg2) = 0.5 * (uGrad(icg1) + uGrad(icg1 + 1));
                vGradSeasurfaceHeight(icg2) = 0.5 * (vGrad(icg1) + vGrad(icg1 + 1));
            }
        }
#pragma omp parallel for
        for (size_t iy = 0; iy < smesh->ny; ++iy) // vertical
        {
            size_t icg1 = (smesh->nx + 1) * iy;
            size_t icg2 = (2 * smesh->nx + 1) * (2 * iy + 1);
            for (size_t ix = 0; ix <= smesh->nx; ++ix, ++icg1, icg2 += 2) {
                uGradSeasurfaceHeight(icg2) = 0.5 * (uGrad(icg1) + uGrad(icg1 + cg1row));
                vGradSeasurfaceHeight(icg2) = 0.5 * (vGrad(icg1) + vGrad(icg1 + cg1row));
            }
        }

        // midpoints
#pragma omp parallel for
        for (size_t iy = 0; iy < smesh->ny; ++iy) // vertical
        {
            size_t icg1 = (smesh->nx + 1) * iy;
            size_t icg2 = (2 * smesh->nx + 1) * (2 * iy + 1) + 1;
            for (size_t ix = 0; ix < smesh->nx; ++ix, ++icg1, icg2 += 2) {
                uGradSeasurfaceHeight(icg2) = 0.25
                    * (uGrad(icg1) + uGrad(icg1 + 1) + uGrad(icg1 + cg1row)
                        + uGrad(icg1 + cg1row + 1));
                vGradSeasurfaceHeight(icg2) = 0.25
                    * (vGrad(icg1) + vGrad(icg1 + 1) + vGrad(icg1 + cg1row)
                        + vGrad(icg1 + cg1row + 1));
            }
        }
    }
}

template <int DGadvection> void CGDynamicsKernel<DGadvection>::prepareIteration(const DataMap& data)
{
    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(*smesh, cgH, data.at(hiceName));
    VectorManipulations::CGAveragePeriodic(*smesh, cgH);
    Interpolations::DG2CG(*smesh, cgA, data.at(ciceName));
    VectorManipulations::CGAveragePeriodic(*smesh, cgA);

    // Reinit the gradient of the sea surface height. Not done by
    // DataMap as seaSurfaceHeight is always dG(0)
    ComputeGradientOfSeaSurfaceHeight(DynamicsKernel<DGadvection, DGstressComp>::seaSurfaceHeight);

    // limit A to [0,1] and H to [0, ...)
    cgA = cgA.cwiseMin(1.0);
    cgA = cgA.cwiseMax(1.e-4);
    cgH = cgH.cwiseMax(1.e-4);
}

template <int CG>
Eigen::Matrix<double, CGDOFS(CG), 1> cgLocal(
    const CGVector<CG>& globalVelocity, int cgi, int cgShift);

template <>
Eigen::Matrix<double, CGDOFS(1), 1> cgLocal(const CGVector<1>& vGlobal, int cgi, int cgShift)
{
    Eigen::Matrix<double, CGDOFS(1), 1> vLocal;
    vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + cgShift), vGlobal(cgi + 1 + cgShift);
    return vLocal;
}

template <>
Eigen::Matrix<double, CGDOFS(2), 1> cgLocal(const CGVector<2>& vGlobal, int cgi, int cgShift)
{
    Eigen::Matrix<double, CGDOFS(2), 1> vLocal;
    vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + 2), vGlobal(cgi + cgShift),
        vGlobal(cgi + 1 + cgShift), vGlobal(cgi + 2 + cgShift), vGlobal(cgi + 2 * cgShift),
        vGlobal(cgi + 1 + 2 * cgShift), vGlobal(cgi + 2 + 2 * cgShift);
    return vLocal;
}

template <int DGadvection> void CGDynamicsKernel<DGadvection>::projectVelocityToStrain()
{
    // !!! must still be converted to the spherical system!!!

    const int cgshift = CGdegree * smesh->nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < smesh->ny; ++row) {
        int dgi = smesh->nx * row; //!< Index of dg vector
        int cgi = CGdegree * cgshift * row; //!< Lower left index of cg vector

        for (size_t col = 0; col < smesh->nx;
             ++col, ++dgi, cgi += CGdegree) { // loop over all elements

            if (smesh->landmask[dgi] == 0) // only on ice
                continue;

            // get the local x/y - velocity coefficients on the element
            Eigen::Matrix<double, CGDOFS(CGdegree), 1> vx_local
                = cgLocal<CGdegree>(u, cgi, cgshift);
            Eigen::Matrix<double, CGDOFS(CGdegree), 1> vy_local
                = cgLocal<CGdegree>(v, cgi, cgshift);

            // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
            // by integrating rhs and inverting with dG(stress) mass matrix
            //
            e11.row(dgi) = pmap->iMgradX[dgi] * vx_local;
            e22.row(dgi) = pmap->iMgradY[dgi] * vy_local;
            e12.row(dgi) = 0.5 * (pmap->iMgradX[dgi] * vy_local + pmap->iMgradY[dgi] * vx_local);

            if (smesh->CoordinateSystem == SPHERICAL) {
                e11.row(dgi) -= pmap->iMM[dgi] * vy_local;
                e12.row(dgi) += 0.5 * pmap->iMM[dgi] * vx_local;
            }
        }
    }
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::addStressTensorCell(
    const size_t eid, const size_t cx, const size_t cy)
{
    Eigen::Vector<Nextsim::FloatType, CGdof> tx = (pmap->divS1[eid] * s11.row(eid).transpose()
        + pmap->divS2[eid] * s12.row(eid).transpose());
    Eigen::Vector<Nextsim::FloatType, CGdof> ty = (pmap->divS1[eid] * s12.row(eid).transpose()
        + pmap->divS2[eid] * s22.row(eid).transpose());

    if (smesh->CoordinateSystem == SPHERICAL) {
        tx += pmap->divM[eid] * s12.row(eid).transpose();
        ty -= pmap->divM[eid] * s11.row(eid).transpose();
    }
    const size_t cgRow = CGdegree * smesh->nx + 1;
    const size_t cg_i
        = CGdegree * cgRow * cy + CGdegree * cx; //!< lower left CG-index in element (cx,cy)

    // Fill the stress divergence values
    for (size_t row = 0; row <= CGdegree; ++row) {
        for (size_t col = 0; col <= CGdegree; ++col) {
            dStressX(cg_i + col + row * cgRow) -= tx(col + (CGdegree + 1) * row);
            dStressY(cg_i + col + row * cgRow) -= ty(col + (CGdegree + 1) * row);
        }
    }
}

template <int DGadvection> void CGDynamicsKernel<DGadvection>::stressDivergence()
{
    // Somewhat meaningless, but it uses the name in the former version of the code
    auto& tx = dStressX;
    auto& ty = dStressY;

    // Zero the stress gradient vectors
#pragma omp parallel for
    for (size_t i = 0; i < tx.rows(); ++i) {
        tx(i) = 0.0;
        ty(i) = 0.0;
    }

    // parallelization in stripes
    for (size_t p = 0; p < 2; ++p) {
#pragma omp parallel for schedule(static)
        for (size_t cy = 0; cy < smesh->ny; ++cy) {
            //!< loop over all cells of the mesh
            if (cy % 2 == p) {
                size_t c = smesh->nx * cy;
                for (size_t cx = 0; cx < smesh->nx; ++cx, ++c)
                    //!< loop over all cells of the mesh
                    if (smesh->landmask[c] == 1) // only on ice!
                        addStressTensorCell(c, cx, cy);
            }
        }
    }
    // set zero on the Dirichlet boundaries
    dirichletZero(tx);
    dirichletZero(ty);
    // add the contributions on the periodic boundaries
    VectorManipulations::CGAveragePeriodic(*smesh, tx);
    VectorManipulations::CGAveragePeriodic(*smesh, ty);
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::dirichletZero(CGVector<CGdegree>& v) const
{
    // the four segments bottom, right, top, left, are each processed in parallel
    for (size_t seg = 0; seg < 4; ++seg) {
#pragma omp parallel for
        for (size_t i = 0; i < smesh->dirichlet[seg].size(); ++i) {

            const size_t eid = smesh->dirichlet[seg][i];
            const size_t ix = eid % smesh->nx; // compute coordinates of element
            const size_t iy = eid / smesh->nx;

            if (seg == 0) // bottom
                for (size_t j = 0; j < CGdegree + 1; ++j)
                    v(iy * CGdegree * (CGdegree * smesh->nx + 1) + CGdegree * ix + j, 0) = 0.0;
            else if (seg == 1) // right
                for (size_t j = 0; j < CGdegree + 1; ++j)
                    v(iy * CGdegree * (CGdegree * smesh->nx + 1) + CGdegree * ix + CGdegree
                            + (CGdegree * smesh->nx + 1) * j,
                        0)
                        = 0.0;
            else if (seg == 2) // top
                for (size_t j = 0; j < CGdegree + 1; ++j)
                    v((iy + 1) * CGdegree * (CGdegree * smesh->nx + 1) + CGdegree * ix + j, 0)
                        = 0.0;
            else if (seg == 3) // left
                for (size_t j = 0; j < CGdegree + 1; ++j)
                    v(iy * CGdegree * (CGdegree * smesh->nx + 1) + CGdegree * ix
                            + (CGdegree * smesh->nx + 1) * j,
                        0)
                        = 0.0;
            else {
                std::cerr << "That should not have happened!" << std::endl;
                abort();
            }
        }
    }
}

template <int DGadvection> void CGDynamicsKernel<DGadvection>::applyBoundaries()
{
    dirichletZero(u);
    dirichletZero(v);
    // TODO Periodic boundary conditions.
}

template <int DGadvection>
std::pair<CGVector<CGdegree>, CGVector<CGdegree>>
CGDynamicsKernel<DGadvection>::getIceOceanStress() const
{
    CGVector<CGdegree> taux, tauy;
    taux.resizeLike(uStress);
    tauy.resizeLike(vStress);

#pragma omp parallel for
    for (int i = 0; i < taux.rows(); ++i) {
        const double uOceanRel = uOcean(i) - uStress(i);
        const double vOceanRel = vOcean(i) - vStress(i);
        const double cPrime = DynamicsKernel<DGadvection, DGstressComp>::getParams().F_ocean
            * std::hypot(uOceanRel, vOceanRel);

        taux(i) = cPrime * (uOceanRel * cosOceanAngle - vOceanRel * sinOceanAngle);
        tauy(i) = cPrime * (vOceanRel * cosOceanAngle + uOceanRel * sinOceanAngle);
    }

    return { taux, tauy };
}

// Instantiate the templates for all (1, 3, 6) degrees of DGadvection
template class CGDynamicsKernel<1>;
template class CGDynamicsKernel<3>;
template class CGDynamicsKernel<6>;

}
