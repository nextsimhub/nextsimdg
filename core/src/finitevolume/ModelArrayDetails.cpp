/*!
 * @file ModelArrayDetails.cpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

// A source file to detail the ModelArray dimensions and types for
// finite volume models, as well as the relationships between them.

// Should be grouped with a consistent ModelArrayDetails.cpp and
// ModelArrayTypedefs.hpp

namespace Nextsim {
std::map<ModelArray::Dimension, ModelArray::DimensionSpec> ModelArray::definedDimensions = {
    { ModelArray::Dimension::X, { "x", 0 } },
    { ModelArray::Dimension::Y, { "y", 0 } },
    { ModelArray::Dimension::Z, { "z", 1 } },
};

ModelArray::TypeDimensions ModelArray::typeDimensions = {
    { ModelArray::Type::H,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::U,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::V,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::Z,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
            ModelArray::Dimension::Z,
        } },
};

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::H, "HField" },
    { ModelArray::Type::U, "UField" },
    { ModelArray::Type::V, "VField" },
    { ModelArray::Type::Z, "ZField" },
};

ModelArray::ModelArray()
    : ModelArray(Type::H)
{
}

bool ModelArray::hasDoF(const Type type) { return false; }

ModelArray::SizeMap::SizeMap()
    : m_sizes({ { Type::H, 0 }, { Type::U, 0 }, { Type::V, 0 }, { Type::Z, 0 } })
{
}

ModelArray::DimensionMap::DimensionMap()
    : m_dimensions({
        { Type::H, { 0 } },
        { Type::U, { 0 } },
        { Type::V, { 0 } },
        { Type::Z, { 0, 1 } },
    })
{
}

const std::vector<ModelArray::Type, ModelArray::Dimension> ModelArray::componentMap = {};

}
