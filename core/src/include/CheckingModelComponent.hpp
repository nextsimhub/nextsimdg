/*!
 * @file CheckingModelComponent.cpp
 *
 * @date 30 May 2025
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
                field.arrayRef->checkLimits(getOceanMask());
            } catch (const std::exception& e) {
                throw std::runtime_error("Check failed for '" + field.name + "': " + e.what());
            }
        }
    }
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
    std::vector<fieldInfo> fieldsToCheck;

    bool boolCheckFields = false;

    static bool& boolCheckAll()
    {
        static bool checkAll = false;
        return checkAll;
    }
};

} /* namespace Nextsim */

#endif // CHECKINGMODELCOMPONENT_HPP
