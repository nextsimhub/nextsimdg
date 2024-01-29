/*!
 * @file CGDynamicsKernel.cpp
 *
 * @date Jan 26, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

/*
 * The implementation of DynamicsKernel which uses continuous Galerkin (CG) numerics.
 */

#include "include/DynamicsKernel.hpp"

// Import this from the build system *somehow*
static const int CGdegree = 2;

namespace Nextsim {

template <int DGadvection, int DGstress>
DynamicsKernel<DGadvection, DGstress>::DynamicsKernel() = default;

template <int DGadvection, int DGstress>
DynamicsKernel<DGadvection, DGstress>::~DynamicsKernel() = default;

template <int DGadvection, int DGstress>
void DynamicsKernel<DGadvection, DGstress>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    if (isSpherical)
        throw std::runtime_error("DG dynamics do not yet handle spherical coordinates.");
        // TODO handle spherical coordinates

    //! Define the spatial mesh
    smesh = new ParametricMesh(Nextsim::CARTESIAN);

    smesh->coordinatesFromModelArray(coords);
    smesh->landmaskFromModelArray(mask);
    smesh->dirichletFromMask();
    // TODO: handle periodic and open edges
    for (ParametricMesh::Edge edge : ParametricMesh::edges) {
        smesh->dirichletFromEdge(edge);
    }

    //! Initialize transport
    dgtransport = new Nextsim::DGTransport<DGadvection>(*smesh);
    dgtransport->settimesteppingscheme("rk2");

    //! Initialize stress transport
    stresstransport = new Nextsim::DGTransport<CG2DGSTRESS(CGdegree)>(*smesh);
    stresstransport->settimesteppingscheme("rk2");

    //! Initialize momentum
    internals.momentum = new Nextsim::CGParametricMomentum<CGdegree>(*smesh);

    // resize CG and DG vectors
    hice.resize_by_mesh(*smesh);
    cice.resize_by_mesh(*smesh);

    internals.u.resize_by_mesh(*smesh);
    internals.v.resize_by_mesh(*smesh);
}

template <int DGdegree>
class DynamicsInternals {

    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    Nextsim::CGParametricMomentum<CGdegree>* momentum;

    friend DynamicsKernel<DGdegree, CG2DGSTRESS(CGdegree)>;
};
}
