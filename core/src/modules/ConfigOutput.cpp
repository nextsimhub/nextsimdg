/*!
 * @file ConfigOutput.cpp
 *
 * @date Aug 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfigOutput.hpp"
#include "include/Logged.hpp"
#include "include/StructureFactory.hpp"

#include <cmath>
#include <sstream>

namespace Nextsim {

const std::string ConfigOutput::all = "ALL";
const std::string ConfigOutput::defaultLastOutput = "0-01-01T00:00:00Z";

template <>
const std::map<int, std::string> Configured<ConfigOutput>::keyMap = {
    { ConfigOutput::PERIOD_KEY, "ConfigOutput.period" },
    { ConfigOutput::START_KEY, "ConfigOutput.start" },
    { ConfigOutput::FIELDNAMES_KEY, "ConfigOutput.field_names" },
};

ConfigOutput::ConfigOutput()
    : m_filePrefix()
    , outputPeriod()
    , firstOutput(true)
    , everyTS(false)
    , outputAllTheFields(false)
    , lastOutput(defaultLastOutput)
    , fieldsForOutput()
    , currentFileName()
    {
    }

void ConfigOutput::configure()
{
    std::string periodString = Configured::getConfiguration(keyMap.at(PERIOD_KEY), std::string(""));
    if (periodString.empty()) {
        everyTS = true;
    } else {
        outputPeriod.parse(periodString);
    }
    std::string startString = Configured::getConfiguration(keyMap.at(START_KEY), std::string(""));
    if (startString.empty()) {
        // If you start the model before 1st January year 0, tough.
        lastOutput.parse(defaultLastOutput);
    } else {
        lastOutput.parse(startString);
        if (!everyTS) {
            lastOutput -= outputPeriod;
        }
    }

    std::string outputFields
        = Configured::getConfiguration(keyMap.at(FIELDNAMES_KEY), std::string(""));
    if (outputFields == all || outputFields.empty()) { // Output *all* the fields?
        outputAllTheFields = true; // Output all the fields!
    } else {
        std::istringstream fieldStream;
        fieldStream.str(outputFields);
        for (std::string line; std::getline(fieldStream, line, ',');) {
            fieldsForOutput.insert(line);
        }
    }
}

void ConfigOutput::outputState(const ModelState& fullState, const ModelMetadata& meta)
{
    if (currentFileName == "") {
        std::stringstream startStream;
        startStream << meta.time();
        currentFileName = m_filePrefix + ".nc";
    }

    ModelState state;
    const ModelState* pState;
    if (outputAllTheFields) {
        pState = &fullState;
    } else {
        // Filter only the given fields to the output state
        for (auto fieldEntry : fullState.data) {
            if (fieldsForOutput.count(fieldEntry.first) > 0) {
                state.data[fieldEntry.first] = fieldEntry.second;
            }
        }
        pState = &state;
    }

    /*
     * Produce output either:
     *    • on every timestep after the start time initially stored in lastOutput
     *    • whenever the current time is an integer number of time periods from the
     *      last output time.
     */
    if ((everyTS && meta.time() >= lastOutput) || (std::fmod((meta.time() - lastOutput).seconds(), outputPeriod.seconds()) == 0.)) {
        Logged::info("ConfigOutput: Outputting " + std::to_string(pState->data.size()) + " fields to "
                + currentFileName + " at " + meta.time().format() + "\n");
        StructureFactory::fileFromState(*pState, meta, currentFileName, false);
        lastOutput = meta.time();
    }
}
} /* namespace Nextsim */
