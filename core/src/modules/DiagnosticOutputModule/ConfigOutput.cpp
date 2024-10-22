/*!
 * @file ConfigOutput.cpp
 *
 * @date 2 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfigOutput.hpp"
#include "include/FileCallbackCloser.hpp"
#include "include/Logged.hpp"
#include "include/StructureFactory.hpp"

#include <cmath>
#include <regex>
#include <sstream>

namespace Nextsim {

const std::string ConfigOutput::all = "ALL";
const std::string ConfigOutput::defaultLastOutput = "0-01-01T00:00:00Z";

static const std::regex ncSuffix(".nc$");

static const std::string pfx = "ConfigOutput";
static const std::string periodKey = pfx + ".period";
static const std::string startKey = pfx + ".start";
static const std::string fieldNamesKey = pfx + ".field_names";
static const std::string fileNameKey = pfx + ".filename";
static const std::string filePeriodKey = pfx + ".file_period";

// Access the model.start key. There's no clean way of getting this from Model, I think.
static const std::string modelStartKey = "model.start";

template <>
const std::map<int, std::string> Configured<ConfigOutput>::keyMap = {
    { ConfigOutput::PERIOD_KEY, periodKey },
    { ConfigOutput::START_KEY, startKey },
    { ConfigOutput::FIELDNAMES_KEY, fieldNamesKey },
    { ConfigOutput::FILENAME_KEY, fileNameKey },
    { ConfigOutput::FILEPERIOD_KEY, filePeriodKey },
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

ConfigurationHelp::HelpMap& ConfigOutput::getHelpText(HelpMap& map, bool getAll)
{
    map[pfx] = {
        { periodKey, ConfigType::STRING, {}, "", "", "Time between samples of the output data." },
        { startKey, ConfigType::STRING, {}, "model.start", "",
            "Date at which to start outputting data." },
        { fieldNamesKey, ConfigType::STRING, {}, "ALL", "",
            "Comma separated, space free list of fields to be output. "
            "The special value \""
                + all + "\" will output all available fields." },
        { fileNameKey, ConfigType::STRING, {}, "", "",
            "Filename pattern for the output diagnostic files. Date and time elements can be "
            "included as in std::put_time()." },
        { filePeriodKey, ConfigType::STRING, {}, "", "",
            "The period with which diagnostic files are created." },
    };
    return map;
}

ConfigurationHelp::HelpMap& ConfigOutput::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    return map;
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
        startString = Configured::getConfiguration(modelStartKey, std::string(""));
        if (startString.empty())
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

    std::string rawFileName = Configured::getConfiguration(fileNameKey, m_filePrefix);
    // Strip any ".nc" suffix from the end of the configured value.
    std::smatch match;
    std::regex_search(rawFileName, match, ncSuffix);
    m_filePrefix = match.empty() ? rawFileName : match.prefix();

    // The default string is the number of seconds in 10000 years of 365 days
    std::string newFilePeriodStr
        = Configured::getConfiguration(keyMap.at(FILEPERIOD_KEY), std::string("315360000000"));
    fileChangePeriod = Duration(newFilePeriodStr);
    lastFileChange = lastOutput;
}

void ConfigOutput::setModelStart(const TimePoint& modelStart)
{
    // Set the lastOutput time to the model start if the default value has not yet been replaced.
    if (lastOutput == TimePoint(defaultLastOutput)) {
        lastOutput = modelStart;
    }
}

void ConfigOutput::outputState(const ModelMetadata& meta)
{
    const TimePoint& time = meta.time();
    if (currentFileName == "" || (lastFileChange + fileChangePeriod <= time)) {
        std::string newFileName = time.format(m_filePrefix) + ".nc";
        if (newFileName != currentFileName) {
            // TODO: Close the file currentFileName
            FileCallbackCloser::close(currentFileName);
            currentFileName = newFileName;
        }
        lastFileChange = time;
    }

    ModelState state;
    auto storeData = ModelComponent::getStore().getAllData();
    if (outputAllTheFields) {
        // If the internal to external name lookup table is still empty, fill it
        if (reverseExternalNames.empty()) {
            for (auto entry : externalNames) {
                // Add the reverse lookup between external and internal names, if one has not been
                // added
                if (!reverseExternalNames.count(entry.second)) {
                    reverseExternalNames[entry.second] = entry.first;
                }
            }
        }

        // Output every entry in storeData, as either its external name if
        // defined, or as its internal name.
        for (auto entry : storeData) {
            if (entry.second && entry.second->trueSize()) {
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
            if (externalNames.count(fieldExtName) && storeData.count(externalNames.at(fieldExtName))
                && storeData.at(externalNames.at(fieldExtName))) {
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
    Duration timeSinceOutput = meta.time() - lastOutput;
    if (timeSinceOutput.seconds() > 0
        && (everyTS || std::fmod(timeSinceOutput.seconds(), outputPeriod.seconds()) == 0.)) {
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
