/*!
 * @file CheckingModelComponent.cpp
 *
 * @date 26 Feb 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/CheckingModelComponent.hpp"

#include <include/PhysicalBounds.hpp>

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
        const ModelArray::DataType masked
            = getOceanMask().data().select(field.arrayRef->data(), fillValue);

        if (masked.isNaN().any()) {
            throw std::runtime_error(
                field.name + " contains a NaN. Error at time step " + tst.start.format());
        }
        if (masked.minCoeff() < field.bounds.first) {
            int i;
            throw std::runtime_error(field.name + " is out of bounds. "
                + std::to_string(masked.col(0).minCoeff(&i)) + " < "
                + std::to_string(field.bounds.first) + ". Error at index " + std::to_string(i)
                + " and time step " + tst.start.format());
        }
        if (masked.maxCoeff() > field.bounds.second) {
            int i;
            throw std::runtime_error(field.name + " is out of bounds. "
                + std::to_string(masked.col(0).maxCoeff(&i)) + " > "
                + std::to_string(field.bounds.second) + ". Error at index " + std::to_string(i)
                + " and time step " + tst.start.format());
        }
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

} /* namespace Nextsim */
