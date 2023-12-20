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
    { ModelArray::Dimension::XVERTEX, { "xvertex", 1 } }, // defined as x + 1
    { ModelArray::Dimension::YVERTEX, { "yvertex", 1 } }, // defined as y + 1
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
        { ModelArray::Type::VERTEX,
            {
                ModelArray::Dimension::XVERTEX,
                ModelArray::Dimension::YVERTEX,
            } },
};

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::H, "HField" },
    { ModelArray::Type::U, "UField" },
    { ModelArray::Type::V, "VField" },
    { ModelArray::Type::Z, "ZField" },
    { ModelArray::Type::VERTEX, "VertexField" },
};

ModelArray::ModelArray()
    : ModelArray(Type::H)
{
}

bool ModelArray::hasDoF(const Type type) { return type == Type::VERTEX; }

ModelArray::SizeMap::SizeMap()
    : m_sizes({ { Type::H, 0 }, { Type::U, 0 }, { Type::V, 0 }, { Type::Z, 0 },  { Type::VERTEX, 1 }, })
{
}

ModelArray::DimensionMap::DimensionMap()
    : m_dimensions({
        { Type::H, { 0 } },
        { Type::U, { 0 } },
        { Type::V, { 0 } },
        { Type::Z, { 0, 1 } },
        { Type::VERTEX, { 1, 1 } },
    })
{
}

const size_t ModelArray::nCoords = 2;

const std::map<ModelArray::Type, ModelArray::Dimension> ModelArray::componentMap = {
        { Type::VERTEX, Dimension::NCOORDS },
};

}
