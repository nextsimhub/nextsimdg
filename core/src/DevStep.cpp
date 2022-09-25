/*!
 * @file DevStep.cpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevStep.hpp"

#include "include/ConfiguredModule.hpp"
#include "include/DiagnosticOutputModule.hpp"

namespace Nextsim {

void DevStep::init()
{
    IDiagnosticOutput& ido = Module::getImplementation<IDiagnosticOutput>();
    // FIXMEL MPI
    ido.setFilenamePrefix("diagnostic");
    tryConfigure(ido);
}

void DevStep::iterate(const TimestepTime& tst)
{
    pData->update(tst);
    // The state of the model has now advanced by one timestep, so update the
    // model metadata timestamp.
    mData->incrementTime(tst.step);
    // XIOS wants all the fields, every timestep, so I guess that's what everyone gets
    ModelState overallState = pData->getStateRecursive(true);
    overallState.merge(ConfiguredModule::getAllModuleConfigurations());
    Module::getImplementation<IDiagnosticOutput>().outputState(overallState, *mData);
}

} /* namespace Nextsim */
