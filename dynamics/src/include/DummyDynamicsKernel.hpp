/*!
 * @file DummyDynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef DUMMYDYNAMICSKERNEL_HPP
#define DUMMYDYNAMICSKERNEL_HPP

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "DGTransport.hpp"

#include "cgVector.hpp"
#include "dgVector.hpp"

#include "CGModelArray.hpp"
#include "DGModelArray.hpp"
#include "include/gridNames.hpp"
#include "include/Time.hpp"


#include <string>
#include <unordered_map>

//Nextsim::COORDINATES CoordinateSystem = Nextsim::CARTESIAN;


namespace Nextsim {

template <int CGdegree, int DGadvection> class DummyDynamicsKernel {
public:
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
            CGModelArray::ma2cg(data, u);
        } else if (name == "v") {
            CGModelArray::ma2cg(data, v);
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

    void update(const TimestepTime& tst) {
        
        const size_t NX = 128;
        const double dt_adv = 120.; //!< Time step of advection problem
        
        //! Define the spatial mesh
        //tmp_create_mesh("tmp-benchmark.smesh", NX, 0.0);
        Nextsim::ParametricMesh smesh(Nextsim::CARTESIAN);
        smesh.readmesh("tmp-benchmark.smesh"); // file temporary committed


        //! Initialize transport
        Nextsim::DGTransport<DGadvection> dgtransport(smesh);
        dgtransport.settimesteppingscheme("rk2");

        //! interpolates CG velocity to DG and reinits normal velocity
        dgtransport.prepareAdvection(u, v);

        //! Perform transport step
        dgtransport.step(dt_adv, cice);
	    dgtransport.step(dt_adv, hice);


    };

private:
    DGVector<DGadvection> hice;
    DGVector<DGadvection> cice;
    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    std::unordered_map<std::string, DGVector<DGadvection>> advectedFields;

    // A map from field name to the type of
    const static std::unordered_map<std::string, ModelArray::Type> fieldType;
};

} /* namespace Nextsim */

#endif /* DUMMYDYNAMICSKERNEL_HPP */
