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
    Module::getImplementation<IDiagnosticOutput>().setFilenamePrefix("diagnostic");
}

void DevStep::iterate(const TimestepTime& tst)
{
    // Get the state of the ocean and atmosphere
    pData->updateAtmosOceanState(tst);
    // Update the model state
    pData->stepPrognosticData(tst);
    // Perform any export of model data not captured in the diagnostic output
    pData->exportData(tst);

    // XIOS wants all the fields, every timestep, so I guess that's what everyone gets
    ModelState overallState = pData->getModelState();
    Module::getImplementation<IDiagnosticOutput>().outputState(overallState, pData->getMetadata());
}

} /* namespace Nextsim */
