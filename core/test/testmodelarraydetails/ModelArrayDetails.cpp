/*!
 * @file ModelArrayDetails.cpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

namespace Nextsim {
std::map<ModelArray::Dimension, ModelArray::DimensionSpec> ModelArray::definedDimensions = {
    { ModelArray::Dimension::X, { "x", 1 } },
    { ModelArray::Dimension::Y, { "y", 1 } },
    { ModelArray::Dimension::Z, { "z", 1 } },
    { ModelArray::Dimension::U, { "u", 1 } },
};

ModelArray::TypeDimensions ModelArray::typeDimensions = {
    { ModelArray::Type::ONED,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
        } },
    { ModelArray::Type::TWOD,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
        } },
    { ModelArray::Type::DOSD,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::Z),
            ModelArray::definedDimensions.at(ModelArray::Dimension::U),
        } },
    { ModelArray::Type::FOURD,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Z),
            ModelArray::definedDimensions.at(ModelArray::Dimension::U),
        } },
};

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::ONED, "OneDField" },
    { ModelArray::Type::TWOD, "TwoDField" },
    { ModelArray::Type::DOSD, "TwoDField" },
    { ModelArray::Type::FOURD, "FourDField" },
};

ModelArray::ModelArray()
    : ModelArray(Type::ONED)
{
}

bool ModelArray::hasDoF(const Type type) { return false; }

ModelArray::SizeMap::SizeMap()
    : m_sizes({
        { Type::ONED, 1 },
        { Type::TWOD, 1 },
        { Type::DOSD, 1 },
        { Type::FOURD, 1 },
    })
{
}

ModelArray::DimensionMap::DimensionMap()
    : m_dimensions({
        { Type::ONED, { 1 } },
        { Type::TWOD, { 1, 1 } },
        { Type::DOSD, { 1, 1 } },
        { Type::FOURD, { 1, 1, 1, 1 } },
    })
{
}

}
