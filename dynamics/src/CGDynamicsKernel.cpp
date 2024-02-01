/*!
 * @file CGDynamicsKernel.cpp
 *
 * @date Jan 26, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

/*
 * The implementation of DynamicsKernel which uses continuous Galerkin (CG) numerics.
 */

#include "include/CGDynamicsKernel.hpp"
#include "include/ModelArray.hpp"

#include "include/cgVector.hpp"
#include "include/Interpolations.hpp"
#include "include/ParametricMap.hpp"
#include "include/VectorManipulations.hpp"

namespace Nextsim {

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    DynamicsKernel<DGadvection, DGstressDegree>::initialise(coords, isSpherical, mask);

    //! Initialize momentum
    momentum = new Nextsim::CGParametricMomentum<CGdegree>(*smesh);

    u.resize_by_mesh(*smesh);
    v.resize_by_mesh(*smesh);
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
        Nextsim::Interpolations::DG2CG(*smesh, momentum->GetAtmx(), utmp);
    } else if (name == vWindName) {
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, momentum->GetAtmy(), vtmp);
    } else if (name == uOceanName) {
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, momentum->GetOceanx(), utmp);
    } else if (name == vOceanName) {
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, momentum->GetOceany(), vtmp);
    } else {
        DynamicsKernel<DGadvection, DGstressDegree>::setData(name, data);
    }
}

template <int DGadvection>
ModelArray CGDynamicsKernel<DGadvection>::getDG0Data(const std::string& name)
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
    } else {
        return DynamicsKernel<DGadvection, DGstressDegree>::getDG0Data(name);
    }
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::prepareAdvection()
{
    dgtransport->prepareTransport(u, v);
    stresstransport->prepareTransport(u, v);
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::prepareIteration(const DataMap& data)
{
    // interpolate ice height and concentration to local cg Variables
    Interpolations::DG2CG(smesh, cgA, data.at(hiceName));
    VectorManipulations::CGAveragePeriodic(smesh, cgA);
    Interpolations::DG2CG(smesh, cgH, data.at(ciceName));
    VectorManipulations::CGAveragePeriodic(smesh, cgH);

    // limit A to [0,1] and H to [0, ...)
    cgA = cgA.cwiseMin(1.0);
    cgA = cgA.cwiseMax(1.e-4);
    cgH = cgH.cwiseMax(1.e-4);
}

template <int CG>
Eigen::Matrix<double, CGDOFS(CG), 1> cgLocal(const CGVector<CG>& globalVelocity, int cgi, int cgShift);

template Eigen::Matrix<double, CGDOFS(1), 1> cgLocal(const CGVector<1>& vGlobal, int cgi, int cgShift)
{
    Eigen::Matrix<double, CGDOFS(1), 1> vLocal;
    vLocal << vGlobal(cgi), vGlobal(cgi + 1),
            vGlobal(cgi + cgShift), vGlobal(cgi + 1 + cgShift);
    return vLocal;
}

template Eigen::Matrix<double, CGDOFS(2), 1> cgLocal(const CGVector<2>& vGlobal, int cgi, int cgShift)
{
    Eigen::Matrix<double, CGDOFS(2), 1> vLocal;
    vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + 2),
            vGlobal(cgi + cgShift), vGlobal(cgi + 1 + cgShift), vGlobal(cgi + 2 + cgShift),
            vGlobal(cgi + 2 * cgShift), vGlobal(cgi + 1 + 2 * cgShift), vGlobal(cgi + 2 + 2 * cgShift);
    return vLocal;
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::projectVelocityToStrain()
{
    // !!! must still be converted to the spherical system!!!

    const int cgshift = CGdegree * smesh->nx + 1; //!< Index shift for each row

    // parallelize over the rows
#pragma omp parallel for
    for (size_t row = 0; row < smesh->ny; ++row) {
      int dgi = smesh->nx * row; //!< Index of dg vector
      int cgi = CGdegree * cgshift * row; //!< Lower left index of cg vector

      for (size_t col = 0; col < smesh->nx; ++col, ++dgi, cgi += CGdegree) { // loop over all elements

    if (smesh->landmask[dgi]==0) // only on ice
      continue;

    // get the local x/y - velocity coefficients on the element
    Eigen::Matrix<double, CGDOFS(CGdegree), 1> vx_local = cgLocal<CGdegree>(u, cgi, cgshift);
    Eigen::Matrix<double, CGDOFS(CGdegree), 1> vy_local = cgLocal<CGdegree>(v, cgi, cgshift);

    // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
    // by integrating rhs and inverting with dG(stress) mass matrix
    //
    e11.row(dgi) = pmap.iMgradX[dgi] * vx_local;
    e22.row(dgi) = pmap.iMgradY[dgi] * vy_local;
    e12.row(dgi) = 0.5 * (pmap.iMgradX[dgi] * vy_local + pmap.iMgradY[dgi] * vx_local);

    if (smesh->CoordinateSystem == SPHERICAL)
      {
        e11.row(dgi) -= pmap.iMM[dgi] * vy_local;
        e12.row(dgi) += 0.5 * pmap.iMM[dgi] * vx_local;
      }
      }
    }

}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::addStressTensorCell(const double scale, const size_t eid, const size_t cx, const size_t cy)
{
    Eigen::Vector<Nextsim::FloatType, CGdof> tx = scale * (pmap.divS1[eid] * S11.row(eid).transpose() + pmap.divS2[eid] * S12.row(eid).transpose());
    Eigen::Vector<Nextsim::FloatType, CGdof> ty = scale * (pmap.divS1[eid] * S12.row(eid).transpose() + pmap.divS2[eid] * S22.row(eid).transpose());

    if (smesh->CoordinateSystem == SPHERICAL) {
        tx += scale * pmap.divM[eid] * s12.row(eid).transpose();
        ty -= scale * pmap.divM[eid] * s11.row(eid).transpose();
    }
    const size_t cgRow = CGdegree * smesh.nx + 1;
    const size_t cg_i = CGdegree * cgRow * cy + CGdegree * cx; //!< lower left CG-index in element (cx,cy)

    // Fill the stress divergence values
    for (size_t row = 0; row <= CGdegree; ++row) {
        for (size_t col = 0; col <= CGdegree; ++col) {
            dStressX(cg_i + col + row * cgRow) -= tx(col + (CGdegree +1) * row);
            dStressY(cg_i + col + row * cgRow) -= ty(col + (CGdegree +1) * row);
        }
    }
}

template <int DGadvection>
void CGDynamicsKernel<DGadvection>::calculateStressDivergence(const double scale)
{
    // Somewhat meaningless, but it uses the name in the former version of the code
    auto& tx = dStressX;
    auto& ty = dStressY;

    // Zero the stress gradient vectors
#pragma omp parallel for
    for (size_t i = 0; i < tx.rows(); ++i) {
        tx(i)=0.0;
        ty(i)=0.0;
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
                        addStressTensorCell(scale, c, cx, cy);
            }
        }
    }
    // set zero on the Dirichlet boundaries
    DirichletZero(tx);
    DirichletZero(ty);
    // add the contributions on the periodic boundaries
    VectorManipulations::CGAveragePeriodic(smesh, tx);
    VectorManipulations::CGAveragePeriodic(smesh, ty);
}

// Instantiate the templates for all (1, 2) degrees of DGadvection

template class CGDynamicsKernel<1>;
template void CGDynamicsKernel<1>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask);
template void CGDynamicsKernel<1>::prepareAdvection();
template void CGDynamicsKernel<1>::projectVelocityToStrain();
template class DynamicsInternals<1>;

template class CGDynamicsKernel<2>;
template void CGDynamicsKernel<2>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask);
template void CGDynamicsKernel<2>::prepareAdvection();
template void CGDynamicsKernel<2>::projectVelocityToStrain();
template class DynamicsInternals<2>;

}
