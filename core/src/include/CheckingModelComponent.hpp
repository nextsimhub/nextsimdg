/*!
 * @file CheckingModelComponent.hpp
 *
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CHECKINGMODELCOMPONENT_HPP
#define CHECKINGMODELCOMPONENT_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

class CheckingModelComponent : public ModelComponent {

protected:
    CheckingModelComponent() = default;

    /*!
     * @brief Check fields listed in fieldsToCheck. Throw a runtime_error if values are outside
     * bounds.
     */
    void checkFields() const;

    /*!
     * @brief Add all the fields listed to the vector fieldsToCheck so they will be checked by
     * checkFields()
     *
     * @param fieldsToAdd An std::map of string-s and ModelArrayAccessors-s with field names and
     * references to check.
     */
    void addChecks(const std::map<const std::string, ModelArrayAccessorBase<RO>>& fieldsToAdd);

    //! @param checkFast Set to true to do a quick check of the main prognostics in PrognosticData
    bool checkFast = false;

    /*!
     * @param checkAll Set to true to do an exhaustive check of all the variable registered in the
     * central registry and any additional variable registered by the modules. */
    static bool& checkAll()
    {
        static bool checkAll = false;
        return checkAll;
    }

private:
    class FieldInfo {
    public:
        FieldInfo(std::string n, const ModelArrayAccessorBase<RO>& a)
            : name(std::move(n))
            , arrayRef(a)
        {
        }

        const std::string name;
        ModelArrayAccessorBase<RO> arrayRef;
    };

    std::vector<FieldInfo> fieldsToCheck {};
};

} /* namespace Nextsim */

#endif // CHECKINGMODELCOMPONENT_HPP
