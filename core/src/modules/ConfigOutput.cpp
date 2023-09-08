/*!
 * @file ConfigOutput.cpp
 *
 * @date 7 Sep 2023
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
            if (externalNames.count(fieldName)) {
                internalFieldsForOutput.insert(externalNames.at(fieldName));
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
    auto storeData = ModelComponent::getStore().getAllData();
    if (outputAllTheFields) {
        // If the internal to external name lookup table is still empty, fill it
        if (reverseExternalNames.empty()) {
            for (auto entry : externalNames) {
                // Add the reverse lookup between external and internal names, if one has not been added
                if (!reverseExternalNames.count(entry.second)) {
                    reverseExternalNames[entry.second] = entry.first;
                }
            }
        }

        // Output every entry in storeData, as either its external name if
        // defined, or as its internal name.
        for (auto entry : storeData) {
            if (entry.second) {
                if (reverseExternalNames.count(entry.first)) {
                    state.data[reverseExternalNames.at(entry.first)] = *entry.second;
                } else {
                    state.data[entry.first] = *entry.second;
                }
            }
        }
    } else {
        // Filter only the given fields to the output state
        for (const auto& fieldExtName : fieldsForOutput) {
            if (externalNames.count(fieldExtName) && storeData.count(externalNames.at(fieldExtName)) && storeData.at(externalNames.at(fieldExtName))) {
                    state.data[fieldExtName] = *storeData.at(externalNames.at(fieldExtName));
            }
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
