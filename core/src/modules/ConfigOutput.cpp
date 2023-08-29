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
    : IDiagnosticOutput()
    , m_filePrefix()
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
        // Sort through the list of fields and create lists of Shared- or ProtectedArrays that
        // correspond to the fields.
        for (const std::string& fieldName : fieldsForOutput) {
            if (sharedExternalNames.count(fieldName)) {
                sharedArraysForOutput.insert(
                    sharedArrayNames.at(sharedExternalNames.at(fieldName)));
            } else if (protectedExternalNames.count(fieldName)) {
                protectedArraysForOutput.insert(
                    protectedArrayNames.at(protectedExternalNames.at(fieldName)));
            } else {
                Logged::warning(
                    "ConfigOutput: No field with the name \"" + fieldName + "\" was found.");
            }
        }
    }
}

void ConfigOutput::outputState(const ModelMetadata& meta)
{
    if (currentFileName == "") {
        std::stringstream startStream;
        startStream << meta.time();
        currentFileName = m_filePrefix + ".nc";
    }

    ModelState state;
    if (outputAllTheFields) {
        for (const auto& entry : protectedArrayNames) {
            ModelArrayConstReference macr
                = getProtectedArray().at(static_cast<size_t>(entry.second));
            if (macr && macr->trueSize() > 0)
                state.data[entry.first] = *macr;
        }
        for (const auto& entry : sharedArrayNames) {
            ModelArrayReference mar = getSharedArray().at(static_cast<size_t>(entry.second));
            if (mar && mar->trueSize() > 0)
                state.data[entry.first] = *mar;
        }
    } else {
        // Filter only the given fields to the output state
        for (const auto& fieldExtName : fieldsForOutput) {
            if (protectedExternalNames.count(fieldExtName)) {
                ModelArrayConstReference macr = getProtectedArray().at(static_cast<size_t>(
                    protectedArrayNames.at(protectedExternalNames.at(fieldExtName))));
                if (macr)
                    state.data[fieldExtName] = *macr;
            } else if (sharedExternalNames.count(fieldExtName)) {
                ModelArrayReference mar = getSharedArray().at(
                    static_cast<size_t>(sharedArrayNames.at(sharedExternalNames.at(fieldExtName))));
                if (mar)
                    state.data[fieldExtName] = *mar;
            } // else do not add any data to the state under that name
        }
    }

    /*
     * Produce output either:
     *    • on every timestep after the start time initially stored in lastOutput
     *    • whenever the current time is an integer number of time periods from the
     *      last output time.
     */
    if ((everyTS && meta.time() >= lastOutput)
        || (std::fmod((meta.time() - lastOutput).seconds(), outputPeriod.seconds()) == 0.)) {
        Logged::info("ConfigOutput: Outputting " + std::to_string(state.data.size()) + " fields to "
            + currentFileName + " at " + meta.time().format() + "\n");
        StructureFactory::fileFromState(state, meta, currentFileName, false);
        lastOutput = meta.time();
    }
}

std::string concatenateFields(const std::set<std::string>& strSet)
{
    std::string outStr = "";
    for (auto& str : strSet) {
        outStr += str + ",";
    }
    return outStr;
}

ModelState ConfigOutput::getStateRecursive(const OutputSpec& os) const
{
    return { {},
        {
            { keyMap.at(PERIOD_KEY), outputPeriod.format() },
            { keyMap.at(START_KEY), lastOutput.format() }, // FIXME Not necessarily the start date!
            { keyMap.at(FIELDNAMES_KEY), concatenateFields(fieldsForOutput) },
        } };
}

} /* namespace Nextsim */
