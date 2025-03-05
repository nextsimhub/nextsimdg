/*!
 * @file CheckingModelComponent.cpp
 *
 * @date 05 Mar 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/CheckingModelComponent.hpp"

namespace Nextsim {

const std::string CheckingModelComponent::all = "ALL";

CheckingModelComponent::CheckingModelComponent()
    : externalNames({
#include "include/ProtectedArrayNames.ipp"
#include "include/SharedArrayNames.ipp"
      }) {};

void CheckingModelComponent::checkFields(const TimestepTime& tst)
{
    for (auto& field : fieldsToCheck) {
        // Only check arrays that are in use (have been set/resized)
        if (field.arrayRef->data().rows() == 0)
            continue;

        // We need a sensible fill value for land points
        const double fillValue = (field.bounds.second + field.bounds.first) * 0.5;

        // We need to mask the data, so that the checks don't stumble over the land mask
        const auto masked = mask(*field.arrayRef, fillValue);

        // Check first for NaNs. The code is different for the bounds check, because Eigen doesn't
        // return an index for NaN-checking.
        if (masked.data().isNaN().any()) {
            throw std::runtime_error(
                field.name + " contains a NaN. Error at time step " + tst.start.format());
        }

        // Now we check the bounds and set the array index (i) and value if we're out of bounds
        size_t i;
        double value;
        if (masked.data().minCoeff() < field.bounds.first) {
            value = masked.data().col(0).minCoeff(&i);
        } else if (masked.data().maxCoeff() > field.bounds.second) {
            value = masked.data().col(0).maxCoeff(&i);
        } else {
            continue;
        }

        // If we haven't continue'd (or thrown an exception) by now we have an error in the field
        const std::vector<size_t> loc = ModelArray::locationFromIndex(masked.getType(), i);
        std::string locStr = "[";
        for (const size_t& l : loc)
            locStr += std::to_string(l) + ",";
        locStr.pop_back();
        locStr.push_back(']');

        throw std::runtime_error(field.name + " is out of bounds. " + std::to_string(value)
            + " not in [" + std::to_string(field.bounds.first) + ","
            + std::to_string(field.bounds.second) + "]. Error at index " + locStr
            + " and time step " + tst.start.format());
    }
}

// Go through the user supplied list of fields and get ModelArray* for each
void CheckingModelComponent::setFieldsToCheck(
    const std::string& listOfFields, const std::string& prefix)
{
    // Populate a map of field name and pointers. Either through getStore.getAllData() or from a
    // user-supplied list
    std::unordered_map<std::string, const ModelArray*> storeData;
    if (listOfFields == all) { // Check *all* the fields
        storeData = getStore().getAllData();
    } else { // Populate storeData with the fields listed
        std::istringstream fieldStream;
        fieldStream.str(listOfFields);
        for (std::string fieldName; std::getline(fieldStream, fieldName, ',');) {
            try {
                storeData.emplace(externalNames.at(fieldName),
                    getStore().getArrayRef(externalNames.at(fieldName)));
            } catch (std::out_of_range&) {
                Logged::warning(
                    prefix + ": No field with the name \"" + fieldName + "\" was found.");
            }
        }
    }

    // Create a reverse look up for externalNames so the user can see the "user-friendly name"
    // in case of a crash
    std::map<std::string, std::string> reverseNames;
    for (const auto& ext : externalNames) {
        reverseNames[ext.second] = ext.first;
    }

    // Populate fieldsToCheck with user supplied fields
    for (const auto& x : storeData) {
        PhysicalBounds bounds;
        fieldsToCheck.emplace_back(
            prefix + ": " + reverseNames.at(x.first), x.second, bounds.getBounds(x.first));
    }
}

ConfigurationHelp::OptionList CheckingModelComponent::getHelpList(const std::string& fieldNamesKey,
    const std::string& fieldNamesDefault, const std::string& checkFieldsKey,
    const bool checkFieldsDefault)
{
    ConfigurationHelp::OptionList options = {
        { fieldNamesKey, ConfigurationHelp::ConfigType::STRING, {}, fieldNamesDefault, "",
            "Comma separated, space free list of fields to be checked if check_fields is true. "
            "The special value \""
                + all + "\" will output all available fields." },
        { checkFieldsKey, ConfigurationHelp::ConfigType::BOOLEAN, { "true", "false" },
            ConfigurationHelp::toString(checkFieldsDefault), "",
            "Switch to check if the fields listed in field_names are within a reasonable "
            "physical range and not NaN." },
    };
    return options;
}

} /* namespace Nextsim */
