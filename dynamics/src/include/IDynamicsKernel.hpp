/*!
 * @file IDynamicsKernel.hpp
 *
 * @date 3 Nov 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICSKERNEL_HPP
#define IDYNAMICSKERNEL_HPP

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

namespace Nextsim {

template<int CGdegree, int DGadvection> class IDynamicsKernel {
public:
    void initialisation(){}

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
    void setData(const std::string& name, const ModelArray& data){}

    /*!
     * @brief Returns an HField ModelArray containing the DG0 finite volume
     * component of the named dynamics field.
     *
     * @param name the name of the requested field.
     *
     */
    ModelArray getDG0Data(const std::string& name){return ModelArray(ModelArray::Type::H);}

    /*!
     * @brief Returns a DG or DGSTRESS ModelArray containing the full DG data for
     * the named dynamics field.
     *
     * @param name the name of the requested field.
     */
    ModelArray getDGData(const std::string& name){return ModelArray(ModelArray::Type::H);}

    void update(const TimestepTime& tst){}
protected:
    Nextsim::DGTransport<DGadvection>* dgtransport;
    Nextsim::CGParametricMomentum<CGdegree>* momentum;

    DGVector<DGadvection> hice;
    DGVector<DGadvection> cice;

private:

};

}

#endif /* IDYNAMICSKERNEL_HPP */
