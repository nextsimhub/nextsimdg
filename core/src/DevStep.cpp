/*!
 * @file DevStep.cpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevStep.hpp"

#include "include/DiagnosticOutputModule.hpp"

namespace Nextsim {

void DevStep::init()
{
    Module::setImplementation<IDiagnosticOutput>("Nextsim::SimpleOutput");
    Module::getImplementation<IDiagnosticOutput>().setFilename("diagnostic.nc");
}

void DevStep::iterate(const TimestepTime& tst)
{
    pData->update(tst);
// TODO: More fine grained control than "all the fields, every timestep"
    ModelState overallState = pData->getStateRecursive(true);
    Module::getImplementation<IDiagnosticOutput>().outputState(overallState, tst);
}

} /* namespace Nextsim */
