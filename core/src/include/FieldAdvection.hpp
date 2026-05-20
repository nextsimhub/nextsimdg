/*!
 * @file FieldAdvection.hpp
 *
 * @date Jun 5, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FIELDADVECTION_HPP
#define FIELDADVECTION_HPP

#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#ifdef USE_KOKKOS
#include "include/ModelArrayStore.hpp" // needed for DeviceViewMA
#endif

namespace Nextsim {

class FieldAdvection {
public:
    /*!
     * Advects a ModelArray field by one timestep.
     *
     * @param field The field to be advected.
     * @param tst The timestep specification.
     *
     * @returns A reference to the updated ModelArray.
     */
    static ModelArray& advectField(ModelArray& field, const TimestepTime& tst,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity());

#ifdef USE_KOKKOS
    static void advectField(const DeviceViewMA& field, const TimestepTime& tst,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity());

/*    template<typename ExecSpace>
    static void advectField(ExecSpace execSpace, const DeviceViewMA& field, const TimestepTime& tst,
        FloatType lowerLimit = -std::numeric_limits<FloatType>::infinity(),
        FloatType upperLimit = std::numeric_limits<FloatType>::infinity())
    {}*/
#endif
};

} /* namespace Nextsim */

#endif /* FIELDADVECTION_HPP */
