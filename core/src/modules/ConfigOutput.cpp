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

// clang-format off
// Using a macro protects against typos
#define SHAREDEL(NAME) { #NAME, SharedArray::NAME }
// Mapping from canonical SharedArray name to the SharedArray enum value
const std::map<std::string, ModelComponent::SharedArray> ConfigOutput::sharedArrayNames = {
        SHAREDEL(H_ICE), // Updated ice thickness, ice average, m
        SHAREDEL(C_ICE), // Updated ice concentration
        SHAREDEL(H_SNOW), // Updated snow depth, ice average, m
        SHAREDEL(T_ICE), // Updated ice temperatures, ˚C
        SHAREDEL(Q_IA), // Ice to atmosphere heat flux W m⁻²
        SHAREDEL(Q_IC), // Ice conduction heat flux W m⁻²
        SHAREDEL(Q_IO), // Ice to ocean heat flux W m⁻²
        SHAREDEL(Q_OW), // Open water heat flux W m⁻²
        SHAREDEL(DQIA_DT), // Derivative of Qᵢₐ w.r.t. ice surface temperature  W m⁻² K⁻¹
        SHAREDEL(HSNOW_MELT), // Thickness of snow that melted, m
        SHAREDEL(SUBLIM), // Upward sublimation rate kg m⁻² s⁻¹
        SHAREDEL(DELTA_HICE), // Change in sea ice thickness, m
        SHAREDEL(DELTA_CICE), // Change in sea ice concentration
        SHAREDEL(NEW_ICE), // Volume of new ice formed [m]
};
#undef SHAREDEL

// Using a macro protects against typos
#define PROTEL(NAME) { #NAME, ProtectedArray::NAME }
// Mapping from canonical ProtectedArray name to the ProtectedArray enum value
const std::map<std::string, ModelComponent::ProtectedArray> ConfigOutput::protectedArrayNames = {
        PROTEL(H_ICE), // Ice thickness, cell average, m
        PROTEL(C_ICE), // Ice concentration
        PROTEL(H_SNOW), // Snow depth, cell average, m
        PROTEL(T_ICE), // Ice temperature, ˚C
        PROTEL(T_AIR), // Air temperature, ˚C
        PROTEL(DEW_2M), // Dew point at 2 m, ˚C
        PROTEL(P_AIR), // sea level air pressure, Pa
        PROTEL(MIXRAT), // water vapour mass mixing ratio
        PROTEL(SW_IN), // incoming shortwave flux, W m⁻²
        PROTEL(LW_IN), // incoming longwave flux, W m⁻²
        PROTEL(MLD), // mixed layer depth, m
        PROTEL(SNOW), // snow fall, kg m⁻² s⁻¹
        PROTEL(SSS), // sea surface salinity, PSU
        PROTEL(SST), // sea surface temperature ˚C
        PROTEL(EVAP_MINUS_PRECIP), // E-P atmospheric freshwater flux, kg s⁻¹ m⁻²
        PROTEL(ML_BULK_CP), // Mixed layer bulk heat capacity J K⁻¹ m⁻²
        PROTEL(TF), // Ocean freezing temperature, ˚C
        PROTEL(WIND_SPEED), // Wind speed, m s⁻¹
        PROTEL(HTRUE_ICE), // Ice thickness, ice average, m
        PROTEL(HTRUE_SNOW), // Snow thickness, ice average, m
        PROTEL(OCEAN_U), // x(east)-ward ocean current, m s⁻¹
        PROTEL(OCEAN_V), // y(north)-ward ocean current, m s⁻¹
};
#undef PROTEL
// clang-format on

/*
 * Using a pair of .ipp files to allow definitions of the externally visible
 * names of the fields to defined outside of an actual source file, even if the
 * definition file has a slightly odd format.
 */
const std::map<std::string, std::string> ConfigOutput::protectedExternalNames = {
#include "include/ProtectedArrayNames.ipp"
};
const std::map<std::string, std::string> ConfigOutput::sharedExternalNames = {
#include "include/SharedArrayNames.ipp"
};

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
