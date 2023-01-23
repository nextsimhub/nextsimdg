/*!
 * @file ConfigOutput.cpp
 *
 * @date Aug 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfigOutput.hpp"
#include "include/Logged.hpp"
#include "include/StructureFactory.hpp"

#include <sstream>

namespace Nextsim {

const std::string ConfigOutput::all = "ALL";

template <>
const std::map<int, std::string> Configured<ConfigOutput>::keyMap = {
    { ConfigOutput::PERIOD_KEY, "ConfigOutput.period" },
    { ConfigOutput::START_KEY, "ConfigOutput.start" },
    { ConfigOutput::FIELDNAMES_KEY, "ConfigOutput.field_names" },
};

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
        lastOutput.parse("0-01-01T00:00:00Z");
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
    std::stringstream startStream;
    startStream << meta.time();
    std::string timeFileName = m_filePrefix + "." + startStream.str() + ".nc";

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

    Logged::info("ConfigOutput: Outputting " + std::to_string(pState->data.size()) + " fields to "
        + timeFileName + "\n");

    if ((everyTS && meta.time() >= lastOutput) || (meta.time() >= lastOutput + outputPeriod)) {
        StructureFactory::fileFromState(*pState, meta, timeFileName);
        lastOutput = meta.time();
    }
}

ModelState ConfigOutput::getStateRecursive(const OutputSpec& os) const
{
    return { {},
        {
            { keyMap.at(PERIOD_KEY), outputPeriod.format() },
            { keyMap.at(START_KEY), lastOutput.format() }, // FIXME Not necessarily the start date!
            { keyMap.at(FIELDNAMES_KEY), fieldsForOutput[0] }, // FIXME not all the fields
        } };
}

} /* namespace Nextsim */
