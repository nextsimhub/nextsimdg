/*!
 * @file DevStep.cpp
 *
 * @date 2 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevStep.hpp"

#include "include/ConfiguredModule.hpp"
#include "include/IDiagnosticOutput.hpp"
#include "include/NextsimModule.hpp"
namespace Nextsim {

DevStep::DevStep()
    : pData(nullptr)
    , mData(nullptr)
    , m_restartPeriod(0)
{
}

void DevStep::init()
{
    IDiagnosticOutput& ido = Module::getImplementation<IDiagnosticOutput>();
    ido.setFilenamePrefix("diagnostic");
    tryConfigure(ido);
}

void DevStep::start(const TimePoint& startTime)
{
    // Set the last output time for the restart files to the model start time
    lastOutput = startTime;
    // Set the model start time for the diagnostic output files
    Module::getImplementation<IDiagnosticOutput>().setModelStart(startTime);
}

void DevStep::iterate(const TimestepTime& tst)
{
    pData->update(tst);
    // The state of the model has now advanced by one timestep, so update the
    // model metadata timestamp.
    mData->incrementTime(tst.step);
    if ((m_restartPeriod.seconds() > 0) && (mData->time() >= lastOutput + m_restartPeriod)) {
        std::string currentFileName = mData->time().format(m_restartFileName);
        pData->writeRestartFile(currentFileName, *mData);
        lastOutput = mData->time();
    }
    // XIOS wants all the fields, every timestep, so I guess that's what everyone gets
    OutputSpec os; // The default OutputSpec is all fields, but only cell average values
    ModelState overallState = pData->getStateRecursive(os);
    overallState.merge(ConfiguredModule::getAllModuleConfigurations());
    Module::getImplementation<IDiagnosticOutput>().outputState(*mData);
}

void DevStep::setRestartDetails(const Duration& restartPeriod, const std::string& fileName)
{
    /*
     * A restart period of zero means zero intermediate restart files.
     */
    m_restartPeriod = restartPeriod;
    m_restartFileName = fileName;
}

} /* namespace Nextsim */
