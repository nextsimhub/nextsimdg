/*!
 * @file CheckingModelComponent.cpp
 *
 * @date 05 Mar 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CHECKINGMODELCOMPONENT_HPP
#define CHECKINGMODELCOMPONENT_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/ModelComponent.hpp"
#include "include/PhysicalBounds.hpp"

namespace Nextsim {

class CheckingModelComponent : public ModelComponent {

public:
    CheckingModelComponent();

    /*!
     * @brief Check fields listed in fieldsToCheck. Throw a runtime_error if values are outside
     * bounds.
     */
    void checkFields(const TimestepTime& tst);
    class fieldInfo {
    public:
        fieldInfo(std::string n, const ModelArray* a, const std::pair<double, double>& b)
            : name(std::move(n))
            , arrayRef(a)
            , bounds(b)
        {
        }

        const std::string name;
        const ModelArray* arrayRef;
        const std::pair<double, double> bounds;
    };
    std::vector<fieldInfo> fieldsToCheck;

    /*!
     * @brief Populate the internal fieldsToCheck vector from a list of TextTag strings and model
     * component name. Use field names from ProtectedArrayNames.ipp and SharedArrayNames.ipp. The
     * bounds are from PhysicalBounds.hpp.
     */
    void setFieldsToCheck(const std::string& listOfFields, const std::string& prefix);

    static const std::string all;

    static ConfigurationHelp::OptionList getHelpList(const std::string& fieldNamesKey,
        const std::string& fieldNamesDefault, const std::string& checkFieldsKey,
        const bool checkFieldsDefault);

private:
    void setFieldsToCheck();
    std::map<std::string, std::string> externalNames;
};

} /* namespace Nextsim */

#endif // CHECKINGMODELCOMPONENT_HPP
