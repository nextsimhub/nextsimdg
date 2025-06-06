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
    static ModelArray& advectField(ModelArray& field, const TimestepTime& tst, double lowerLimit =
            -std::numeric_limits<double>::infinity(), double upperLimit =
            std::numeric_limits<double>::infinity());
};

} /* namespace Nextsim */

#endif /* FIELDADVECTION_HPP */
