/*!
 * @file FieldAdvection.cpp
 *
 * @date Jun 5, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FieldAdvection.hpp"

#include "include/IDynamics.hpp"
#include "include/NextsimModule.hpp"

namespace Nextsim {

ModelArray& FieldAdvection::advectField(
    ModelArray& field, const TimestepTime& tst, FloatType lowerLimit, FloatType upperLimit)
{
    Module::getImplementation<IDynamics>().advectField(
        tst.step.seconds(), field, lowerLimit, upperLimit);
    return field;
}

#ifdef USE_KOKKOS
void FieldAdvection::advectField(
    const DeviceViewMA& field, const TimestepTime& tst, FloatType lowerLimit, FloatType upperLimit)
{
    Module::getImplementation<IDynamics>().advectField(
        tst.step.seconds(), field, lowerLimit, upperLimit);
}
#endif

} /* namespace Nextsim */
