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

#include "Tools.hpp"
#include "cgVector.hpp"
#include "dgLimit.hpp"
#include "dgVector.hpp"
#include "dgVisu.hpp"

#include "DGModelArray.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"
#include "include/gridNames.hpp"

#include <map>
#include <set>
#include <string>
#include <unordered_map>

namespace Nextsim {

// forward define the class holding the potentially non-DG parts
template <int DGdegree> class DynamicsInternals;

template <int DGadvection, int DGstress> class DynamicsKernel {
public:
    typedef std::pair<const std::string, const DGVector<DGadvection>&> DataMapping;
    typedef std::map<typename DataMapping::first_type, typename DataMapping::second_type> DataMap;

    DynamicsKernel() = default;
    virtual ~DynamicsKernel() = default;
    virtual void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask)
    {
        //! Define the spatial mesh
        smesh = new ParametricMesh((isSpherical) ? Nextsim::SPHERICAL : Nextsim::CARTESIAN);

        smesh->coordinatesFromModelArray(coords);
        if (isSpherical)
            smesh->RotatePoleToGreenland();
        smesh->landmaskFromModelArray(mask);
        smesh->dirichletFromMask();
        // TODO: handle periodic and open edges
        for (ParametricMesh::Edge edge : ParametricMesh::edges) {
            smesh->dirichletFromEdge(edge);
        }

        //! Initialize transport
        dgtransport = new Nextsim::DGTransport<DGadvection>(*smesh);
        dgtransport->settimesteppingscheme("rk2");

        // resize DG vectors
        hice.resize_by_mesh(*smesh);
        cice.resize_by_mesh(*smesh);
        hice.zero();
        cice.zero();

        e11.resize_by_mesh(*smesh);
        e12.resize_by_mesh(*smesh);
        e22.resize_by_mesh(*smesh);
        s11.resize_by_mesh(*smesh);
        s12.resize_by_mesh(*smesh);
        s22.resize_by_mesh(*smesh);
        e11.zero();
        e12.zero();
        e22.zero();
        s11.zero();
        s12.zero();
        s22.zero();
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
    virtual void setData(const std::string& name, const ModelArray& data)
    {

        // Special cases: hice, cice, (damage, stress) <- not yet implemented
        if (name == hiceName) {
            DGModelArray::ma2dg(data, hice);
        } else if (name == ciceName) {
            DGModelArray::ma2dg(data, cice);
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
    virtual ModelArray getDG0Data(const std::string& name) const
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
    virtual ModelArray getDGData(const std::string& name) const
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

    virtual void update(const TimestepTime& tst) { ++stepNumber; }

    void advectionAndLimits(const TimestepTime& tst)
    {
        prepareAdvection();

        // Set the advection time step length for all advection calls
        deltaTAdvection = tst.step.seconds();
        //! Perform transport step
        dgtransport->step(deltaTAdvection, cice);
        dgtransport->step(deltaTAdvection, hice);

        //! Gauss-point limiting
        Nextsim::LimitMax(cice, 1.0);
        Nextsim::LimitMin(cice, 0.0);
        Nextsim::LimitMin(hice, 0.0);
    }

    /*!
     * Advects a ModelArray by one advection timestep.
     *
     * @param field Data for the field to be advected, as a ModelArray
     * @param fieldName Name of the field to be advected.
     */
    ModelArray& advectField(ModelArray& maField, const std::string& fieldName)
    {
        // Check that advection has been prepared, if not prepare it
        if (!isAdvectionReady) {
            prepareAdvection();
            // Check that the necessary calculations happened
            if (!isAdvectionReady) {
                throw std::runtime_error("DynamicsKernel: Failed to prepare advection.");
            }
        }

        std::map<std::string, size_t> advectedLevels;
        size_t n2DFields = 1;
        // Handle ModelArrays with more than 2 dimensions
        if (maField.nDimensions() > 2) {
            // Get the total number of higher dimensions. ModelArray::Type::H is assumed to be the
            // size of a 2D field.
            size_t size2D = ModelArray::size(ModelArray::Type::H);
            size_t n2DFields = maField.size() / size2D;
            for (size_t k = 0; k < n2DFields; ++k) {
                std::string field2DName = fieldName + generateSuffix(k, n2DFields);
                advectedLevels[field2DName] = k;
            }
        } else {
            // Add one advected level, with level = 0
            advectedLevels[fieldName] = 0;
        }

        for (auto entry : advectedLevels) {
            const std::string& name = entry.first;
            const size_t k = entry.second;
            // Create the DGVector to hold the higher DG components, if it doesn't already exist.
            if (advectedFields.count(name) < 1) {
                DGVector<DGadvection> newField;
                newField.resize_by_mesh(*smesh);
                newField.zero();
                // copy to the map of advected fields
                advectedFields[name] = newField;
            }

            DGVector<DGadvection>& dgField = advectedFields[name];
            // Set the DG0 component for the selected level (or whole field if 2D)
            DGModelArray::ma2dg(maField, dgField, k);

            // Perform the advection
            dgtransport->step(deltaTAdvection, dgField);

            // Extract the DG0 data to the selected level (or whole field if 2D)
            DGModelArray::dg2ma(dgField, maField, k);
        }
        return maField;
    }
protected:
    Nextsim::DGTransport<DGadvection>* dgtransport;

    DGVector<DGadvection> hice;
    DGVector<DGadvection> cice;

    //! Vectors storing strain and stress components
    DGVector<DGstress> e11, e12, e22;
    DGVector<DGstress> s11, s12, s22;

    size_t nSteps = 100;

    size_t stepNumber = 0;

    double deltaT;
    // Advection timestep, specifically
    double deltaTAdvection;
    // Has prepareAdvection been called with the current u & v fields?
    bool isAdvectionReady;
    Nextsim::ParametricMesh* smesh;

    virtual void updateMomentum(const TimestepTime& tst) = 0;

    // Pass through functions to the common momentum solver class
    /*!
     * Prepares the stress iteration and momentum equation solution.
     *
     * @param data A map from data name to the DGVector containing that data.
     *             The keys must include H (hiceName), A (ciceName) and can
     *             optionally contain D (damageName).
     */
    virtual void prepareIteration(const DataMap& data) = 0;
    /*!
     * Updates the strain values based on the velocities
     */
    virtual void projectVelocityToStrain() = 0;
    /*!
     * Calculates the divergence of the stress tensor.
     */
    virtual void stressDivergence() = 0;
    /*!
     * Apply Dirichlet and periodic boundary conditions.
     */
    virtual void applyBoundaries() = 0;

    /*!
     * Prepares the transport objects to perform the advection step.
     */
    virtual void prepareAdvection() = 0;

private:
    std::unordered_map<std::string, DGVector<DGadvection>> advectedFields;

    // A map from field name to the type of
    const std::unordered_map<std::string, ModelArray::Type> fieldType;

    // Generates a consistent suffix for multilevel fields
    std::string generateSuffix(size_t i, size_t n) {
        if (n < 2) return "";
        return "_" + std::to_string(i) + "OF" + std::to_string(n);
    }

};

}

#endif /* DYNAMICSKERNEL_HPP */
