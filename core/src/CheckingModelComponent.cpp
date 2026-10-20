/*!
 * @file CheckingModelComponent.cpp
 *
 * @author Einar Ólason <einar.olason@nersc.no>
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/CheckingModelComponent.hpp"

namespace Nextsim {

void CheckingModelComponent::checkFields() const
{
    // Do nothing if checks are not enabled
    if (!checkFast && !checkAll())
        return;

    for (const auto& field : fieldsToCheck) {
        if (std::optional<std::string> e = field.arrayRef.getHostRO().checkLimits(oceanMask())) {
            throw std::runtime_error("Check failed for '" + field.name + "': " + *e);
        }
    }
}

void CheckingModelComponent::addChecks(
    const std::map<const std::string, ModelArrayAccessorBase<RO>>& fieldsToAdd)
{
    for (const auto& field : fieldsToAdd)
        fieldsToCheck.emplace_back(field.first, field.second);
}

}
