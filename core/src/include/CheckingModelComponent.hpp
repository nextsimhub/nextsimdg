/*!
 * @file ModelComponent.cpp
 *
 * @date 16 Feb 2025
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CHECKINGMODELCOMPONENT_HPP
#define CHECKINGMODELCOMPONENT_HPP

#include "include/ModelComponent.hpp"

namespace Nextsim {

class CheckingModelComponent : public ModelComponent {

public:
    CheckingModelComponent();

    /*!
     * @brief Check fields listed in fieldsToCheck. Throw a runtime_error if values are outside
     * bounds.
     */
    void checkFields(const TimestepTime& tst);
    std::vector<std::tuple<std::string, const ModelArray*, std::pair<double, double>>>
        fieldsToCheck;

    /*!
     * @brief Populate the internal fieldsToCheck vector from a list of TextTag strings and model
     * component name. Use field names from ProtectedArrayNames.ipp and SharedArrayNames.ipp. The
     * bounds are from PhysicalBounds.hpp.
     */
    void setFieldsToCheck(const std::string& listOfFields, const std::string& prefix);

    static const std::string all;

private:
    int numToCheck;
    int layerToCheck;
    void checkFieldsElement(size_t i, const TimestepTime& tst) const;

    void setFieldsToCheck();

    std::map<std::string, std::string> externalNames;
};

} /* namespace Nextsim */

#endif // CHECKINGMODELCOMPONENT_HPP
