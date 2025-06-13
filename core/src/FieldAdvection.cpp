/*!
 * @file FieldAdvection.cpp
 *
 * @date Jun 5, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FieldAdvection.hpp"

#include "include/NextsimModule.hpp"
#include "include/IDynamics.hpp"

namespace Nextsim {

ModelArray& FieldAdvection::advectField(ModelArray& field, const TimestepTime& tst, double lowerLimit, double upperLimit)
{
    Module::getImplementation<IDynamics>().advectField( tst.step.seconds(), field, lowerLimit, upperLimit);
    return field;
}


} /* namespace Nextsim */
