/*!
 * @file ModelComponent.cpp
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
    // numToCheck must be in scope for checkFieldElement as well
    for (numToCheck = 0; numToCheck < fieldsToCheck.size(); ++numToCheck) {

        const ModelArray* array = fieldsToCheck[numToCheck].arrayRef;

        // Only check arrays that are in use (have been set/resized)
        if (array->data().rows() == 0)
            continue;

        int nLayers;
        if (array->nDimensions() == 3)
            nLayers = array->dimensions()[2];
        else if (array->nDimensions() == 2)
            nLayers = 1;
        else
            throw std::logic_error(
                "ModelComponent::checkFields expected a field with 2 or 3 dimensions.\n");

        // layerToCheck must be in scope for checkFieldElement as well
        for (layerToCheck = 0; layerToCheck < nLayers; ++layerToCheck)
            overElements(std::bind(&CheckingModelComponent::checkFieldsElement, this,
                             std::placeholders::_1, std::placeholders::_2),
                tst);
    }
}

void CheckingModelComponent::checkFieldsElement(size_t i, const TimestepTime& tst) const
{
    const auto& x = fieldsToCheck[numToCheck];
    const double& value = x.arrayRef->zIndexAndLayer(i, layerToCheck);
    if (value < x.bounds.first || value > x.bounds.second || std::isnan(value))
        throw std::runtime_error(x.name + " is out of bounds: " + std::to_string(value)
            + " is not within [" + std::to_string(x.bounds.first) + ","
            + std::to_string(x.bounds.second) + "].\n" + "Error at time step " + tst.start.format()
            + " and index " + std::to_string(i));
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