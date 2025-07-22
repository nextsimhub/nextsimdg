/*!
 * @file CheckingModelComponent.hpp
 *
 * @date 02 Jun 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CHECKINGMODELCOMPONENT_HPP
#define CHECKINGMODELCOMPONENT_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

class CheckingModelComponent : public ModelComponent {

protected:
    CheckingModelComponent() = default;

    /*!
     * @brief Check fields listed in fieldsToCheck. Throw a runtime_error if values are outside
     * bounds.
     */
    void checkFields() const
    {
        // Do nothing if checks are not enabled
        if (!boolCheckFields && !boolCheckAll())
            return;

        for (const auto& field : fieldsToCheck) {

            try {
                field.arrayRef->checkLimits(oceanMask());
            } catch (const std::exception& e) {
                throw std::runtime_error("Check failed for '" + field.name + "': " + e.what());
            }
        }
    }
    /*!
     * @brief Add a single field to the vector fieldsToCheck so they will be checked by
     * checkFields()
     *
     * @param fieldToAdd An std::pair of string and ModelArrayReference with field name and
     * reference to check.
     */
    void addChecks(const std::pair<const std::string, const ModelArray*>& fieldToAdd)
    {
        fieldsToCheck.emplace_back(fieldToAdd.first, fieldToAdd.second);
    }
    /*!
     * @brief Add all the fields listed to the vector fieldsToCheck so they will be checked by
     * checkFields()
     *
     * @param fieldsToAdd An std::map of string-s and ModelArrayReference-s with field names and
     * references to check.
     */
    void addChecks(const std::map<const std::string, const ModelArray*>& fieldsToAdd)
    {
        for (const auto& field : fieldsToAdd)
            addChecks(field);
    }

    bool boolCheckFields = false;

    static bool& boolCheckAll()
    {
        static bool checkAll = false;
        return checkAll;
    }

private:
    class fieldInfo {
    public:
        fieldInfo(std::string n, const ModelArray* a)
            : name(std::move(n))
            , arrayRef(a)
        {
        }

        const std::string name;
        const ModelArray* arrayRef;
    };

    std::vector<fieldInfo> fieldsToCheck {};
};

} /* namespace Nextsim */

#endif // CHECKINGMODELCOMPONENT_HPP
