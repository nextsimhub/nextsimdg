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
    ido.setFilenamePrefix("diagnostic");
    tryConfigure(ido);
}

void DevStep::iterate(const TimestepTime& tst)
{
    pData->update(tst);
    // The state of the model has now advanced by one timestep, so update the
    // model metadata timestamp.
    mData->incrementTime(tst.step);
    // Output the model state
    Module::getImplementation<IDiagnosticOutput>().outputState(*mData);
}

} /* namespace Nextsim */
