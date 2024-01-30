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
#include "include/ModelArray.hpp"

// Import this from the build system *somehow*
static const int CGdegree = 2;
static const int DGstressDegree = CG2DGSTRESS(CGdegree);


namespace Nextsim {

template <int DGadvection, int DGstress>
void DynamicsKernel<DGadvection, DGstress>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    initialiseDG(coords, isSpherical, mask);

    //! Initialize momentum
    internals.momentum = new Nextsim::CGParametricMomentum<CGdegree>(*smesh);

    internals.u.resize_by_mesh(*smesh);
    internals.v.resize_by_mesh(*smesh);
}

template <int DGadvection, int DGstress>
void DynamicsKernel<DGadvection, DGstress>::setVelocityData(const std::string& name, const ModelArray& data)
{
    if (name == uName) {
        // FIXME take into account possibility to restart form CG
        // CGModelArray::ma2cg(data, u);
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, internals.u, utmp);
    } else if (name == vName) {
        // CGModelArray::ma2cg(data, v);
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, internals.v, vtmp);
    } else if (name == uWindName) {
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, internals.momentum->GetAtmx(), utmp);
    } else if (name == vWindName) {
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, internals.momentum->GetAtmy(), vtmp);
    } else if (name == uOceanName) {
        DGVector<DGadvection> utmp(*smesh);
        DGModelArray::ma2dg(data, utmp);
        Nextsim::Interpolations::DG2CG(*smesh, internals.momentum->GetOceanx(), utmp);
    } else if (name == vOceanName) {
        DGVector<DGadvection> vtmp(*smesh);
        DGModelArray::ma2dg(data, vtmp);
        Nextsim::Interpolations::DG2CG(*smesh, internals.momentum->GetOceany(), vtmp);
    }
}

template <int DGadvection, int DGstress>
ModelArray DynamicsKernel<DGadvection, DGstress>::getVelocityDG0Data(const std::string& name)
{
    if (name == uName) {
        ModelArray data(ModelArray::Type::U);
        DGVector<DGadvection> utmp(*smesh);
        Nextsim::Interpolations::CG2DG(*smesh, utmp, internals.momentum->u);
        return DGModelArray::dg2ma(utmp, data);
    } else if (name == vName) {
        ModelArray data(ModelArray::Type::V);
        DGVector<DGadvection> vtmp(*smesh);
        Nextsim::Interpolations::CG2DG(*smesh, vtmp, internals.momentum->v);
        return DGModelArray::dg2ma(vtmp, data);
    } else {
        ModelArray noData(ModelArray::Type::H);
        noData.resize();
        noData = 0;
        return noData;
    }
}

template <int DGadvection, int DGstress>
void DynamicsKernel<DGadvection, DGstress>::prepareAdvection()
{
    dgtransport->prepareTransport(internals.momentum->u, internals.momentum->v);
    stresstransport->prepareTransport(internals.momentum->u, internals.momentum->v);
}

template <int DGdegree>
class DynamicsInternals {

    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    Nextsim::CGParametricMomentum<CGdegree>* momentum;

    friend DynamicsKernel<DGdegree, CG2DGSTRESS(CGdegree)>;
};

// Instantiate the templates for all (1, 2) degrees of DGadvection

template class DynamicsKernel<1, DGstressDegree>;
template void DynamicsKernel<1, DGstressDegree>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask);
template void DynamicsKernel<1, DGstressDegree>::setVelocityData(const std::string& name, const ModelArray& data);
template ModelArray DynamicsKernel<1, DGstressDegree>::getVelocityDG0Data(const std::string& name);
template void DynamicsKernel<1, DGstressDegree>::prepareAdvection();
template class DynamicsInternals<1>;

template class DynamicsKernel<2, DGstressDegree>;
template void DynamicsKernel<2, DGstressDegree>::initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask);
template void DynamicsKernel<2, DGstressDegree>::setVelocityData(const std::string& name, const ModelArray& data);
template ModelArray DynamicsKernel<2, DGstressDegree>::getVelocityDG0Data(const std::string& name);
template void DynamicsKernel<2, DGstressDegree>::prepareAdvection();
template class DynamicsInternals<2>;

}
