/*!
 * @file DynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef DYNAMICSKERNEL_HPP
#define DYNAMICSKERNEL_HPP

#include "DGTransport.hpp"
#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "cgParametricMomentum.hpp"

#include "Tools.hpp"
#include "cgVector.hpp"
#include "dgLimit.hpp"
#include "dgVector.hpp"
#include "dgVisu.hpp"

#include "CGModelArray.hpp"
#include "DGModelArray.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"
#include "include/gridNames.hpp"

#include <string>
#include <unordered_map>

#include <iostream>

namespace Nextsim {

template <int CGdegree, int DGadvection> class DynamicsKernel {
public:
    void initialisation(const ModelArray& coords, bool isSpherical, const ModelArray& mask)
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

        //! Initialize momentum
        momentum = new Nextsim::CGParametricMomentum<CGdegree>(*smesh);

        // resize CG and DG vectors
        hice.resize_by_mesh(*smesh);
        cice.resize_by_mesh(*smesh);

        u.resize_by_mesh(*smesh);
        v.resize_by_mesh(*smesh);
    }

    /*!
     * @brief Sets the data from a provided ModelArray.
     *
     * @details Given a name and a ModelArray, sets the data associated with that
     * name. In some special cases (hice, cice…) this is a special array used in
     * the dynamics calculations. In all other cases, these are added to the
     * container of name data fields to be advected. The provided ModelArray can be
     * of DG or DGSTRESS type, in which case all components of the DGVector are
     * filled, or any other type which only fills the DG0 finite volume element of
     * the dgVector. The behaviour is exactly that of the ma2dg() function defined
     * in the DGModelArray class.
     *
     * @param name The name of the data field to set.
     * @param data The ModelArray containing the data to be set.
     *
     */
    void setData(const std::string& name, const ModelArray& data)
    {

        // Special cases: hice, cice, (damage, stress) <- not yet implemented
        if (name == hiceName) {
            DGModelArray::ma2dg(data, hice);
        } else if (name == ciceName) {
            DGModelArray::ma2dg(data, cice);
        } else if (name == "u") {
            // FIXME take into account possibility to restart form CG
            // CGModelArray::ma2cg(data, u);
            DGVector<DGadvection> utmp(*smesh);
            DGModelArray::ma2dg(data, utmp);
            Nextsim::Interpolations::DG2CG(*smesh, u, utmp);
        } else if (name == "v") {
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
            // All other fields get shoved in a (labelled) bucket
            DGModelArray::ma2dg(data, advectedFields[name]);
        }
    }

    /*!
     * @brief Returns an HField ModelArray containing the DG0 finite volume
     * component of the named dynamics field.
     *
     * @param name the name of the requested field.
     *
     */
    ModelArray getDG0Data(const std::string& name)
    {
        HField data(ModelArray::Type::H);
        if (name == hiceName) {
            return DGModelArray::dg2ma(hice, data);
        } else if (name == ciceName) {
            return DGModelArray::dg2ma(cice, data);
        } else if (name == uName) {
            DGVector<DGadvection> utmp(*smesh);
            Nextsim::Interpolations::CG2DG(*smesh, utmp, momentum->GetVx());
            return DGModelArray::dg2ma(utmp, data);
        } else if (name == vName) {
            DGVector<DGadvection> vtmp(*smesh);
            Nextsim::Interpolations::CG2DG(*smesh, vtmp, momentum->GetVy());
            return DGModelArray::dg2ma(vtmp, data);
        } else {
            // Any other named field must exist
            return DGModelArray::dg2ma(advectedFields.at(name), data);
        }
    }

    /*!
     * @brief Returns a DG or DGSTRESS ModelArray containing the full DG data for
     * the named dynamics field.
     *
     * @param name the name of the requested field.
     */
    ModelArray getDGData(const std::string& name)
    {

        if (name == hiceName) {
            DGField data(ModelArray::Type::DG);
            data.resize();
            return DGModelArray::dg2ma(hice, data);
        } else if (name == ciceName) {
            DGField data(ModelArray::Type::DG);
            data.resize();
            return DGModelArray::dg2ma(cice, data);
        } else {
            ModelArray::Type type = fieldType.at(name);
            ModelArray data(type);
            data.resize();
            return DGModelArray::dg2ma(advectedFields.at(name), data);
        }
    }

    void update(const TimestepTime& tst)
    {

        static int step_number = 0;

        //! Perform transport step
        dgtransport->prepareAdvection(momentum->GetVx(), momentum->GetVy());

        dgtransport->step(tst.step.seconds(), cice);
        dgtransport->step(tst.step.seconds(), hice);

        //! Gauss-point limiting
        Nextsim::LimitMax(cice, 1.0);
        Nextsim::LimitMin(cice, 0.0);
        Nextsim::LimitMin(hice, 0.0);

        momentum->prepareIteration(hice, cice);
        //! Momentum
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {
            momentum->mEVPStep(VP, NT_evp, alpha, beta, tst.step.seconds(), hice, cice);
        }

        step_number++;
    };

private:
    DGVector<DGadvection> hice;
    DGVector<DGadvection> cice;
    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    Nextsim::DGTransport<DGadvection>* dgtransport;
    Nextsim::CGParametricMomentum<CGdegree>* momentum;

    Nextsim::ParametricMesh* smesh;

    //! Rheology-Parameters
    Nextsim::VPParameters VP;
    //! MEVP parameters
    double alpha = 1500.0;
    double beta = 1500.0;
    size_t NT_evp = 100;

    std::unordered_map<std::string, DGVector<DGadvection>> advectedFields;

    // A map from field name to the type of
    const static std::unordered_map<std::string, ModelArray::Type> fieldType;
};

} /* namespace Nextsim */

#endif /* DYNAMICSKERNEL_HPP */
