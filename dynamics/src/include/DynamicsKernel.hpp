/*!
 * @file DynamicsKernel.hpp
 *
 * @date Jan 5, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
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

#include "DGModelArray.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"
#include "include/gridNames.hpp"

#include <string>
#include <map>
#include <unordered_map>

namespace Nextsim {

// forward define the class holding the potentially non-DG parts
template <int DGdegree>
class DynamicsInternals;

template <int DGadvection, int DGstress> class DynamicsKernel {
public:

    typedef std::pair<const std::string, const DGVector<DGadvection>&> DataMapping;
    typedef std::map<DataMapping::first_type, DataMapping::second_type> DataMap;

    DynamicsKernel() = default;
    virtual ~DynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask);
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
    virtual void setData(const std::string& name, const ModelArray& data)
    {

        // Special cases: hice, cice, (damage, stress) <- not yet implemented
        if (name == hiceName) {
            DGModelArray::ma2dg(data, hice);
        } else if (name == ciceName) {
            DGModelArray::ma2dg(data, cice);
        } else if (velocityFields.count(name)) {
            setVelocityData(name, data);
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
        } else if (name == uName || name == vName) {
            return getVelocityDG0Data(name);
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
        static int stepNumber = 0;

        prepareAdvection();

        //! Perform transport step
        dgtransport->step(tst.step.seconds(), cice);
        dgtransport->step(tst.step.seconds(), hice);

        //! Perform transport step for stress
        // stresstransport->step(tst.step.seconds(), momentum->GetS11());
        // stresstransport->step(tst.step.seconds(), momentum->GetS12());
        // stresstransport->step(tst.step.seconds(), momentum->GetS22());

        //! Gauss-point limiting
        Nextsim::LimitMax(cice, 1.0);
        Nextsim::LimitMin(cice, 0.0);
        Nextsim::LimitMin(hice, 0.0);

        updateMomentum(tst);

        ++stepNumber;
    }

protected:
    Nextsim::DGTransport<DGadvection>* dgtransport;
    Nextsim::DGTransport<DGstress>* stresstransport;

    DGVector<DGadvection> hice;
    DGVector<DGadvection> cice;

    //! Vectors storing strain and stress components
    DGVector<DGstress> e11, e12, e22;
    DGVector<DGstress> s11, s12, s22;

    size_t NT_evp = 100;

    virtual void updateMomentum(const TimestepTime& tst) = 0;

    // Pass through functions to the common momentum solver class
    /*!
     * Prepares the stress iteration and momentum equation solution.
     *
     * @param data A map from data name to the DGVector containing that data.
     *             The keys must include H (hiceName), A (ciceName) and can
     *             optionally contain D (damageName).
     */
    virtual void prepareIteration(const DataMap& data);
    /*!
     * Updates the strain values based on the velocities
     */
    void projectVelocityToStrain();
    /*!
     * Calculates the divergence of the stress tensor.
     */
    void calculateStressDivergence(const double scale);
    /*!
     * Apply Dirichlet and periodic boundary conditions.
     */
    void applyBoundaries();

    /*!
     * Prepares the transport objects to perform the advection step.
     */
    void prepareAdvection();

    //! Sets the data for a velocity component
    void setVelocityData(const std::string& velocityFieldName, const ModelArray& data);

    //! Gets the DG0/finite volume data for velocity fields
    ModelArray getVelocityDG0Data(const std::string& name);

private:

    DynamicsInternals<DGadvection>& internals;

    Nextsim::ParametricMesh* smesh;

    std::unordered_map<std::string, DGVector<DGadvection>> advectedFields;

    // A map from field name to the type of
    const std::unordered_map<std::string, ModelArray::Type> fieldType;

    // Fields that are specially treated on set and get
    inline static const std::set<std::string> velocityFields = { uName, vName, uWindName, vWindName, uOceanName, vOceanName };

    // Initialises the common DG parts that do not depend on the velocity numerics implementation
    void initialiseDG(const ModelArray& coords, bool isSpherical, const ModelArray& mask)
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
        stresstransport = new Nextsim::DGTransport<DGstress>(*smesh);
        stresstransport->settimesteppingscheme("rk2");

        // resize DG vectors
        hice.resize_by_mesh(*smesh);
        cice.resize_by_mesh(*smesh);
    }
};

}

#endif /* DYNAMICSKERNEL_HPP */
