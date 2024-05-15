/*!
 * @file ModelArrayDetails.cpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

// A source file to detail the ModelArray dimensions and types for
// Discontinuous Galerkin models, as well as the relationships between them.

// Should be grouped with a consistent ModelArrayDetails.cpp and
// ModelArrayTypedefs.hpp

namespace Nextsim {
// clang-format off
std::map<ModelArray::Dimension, ModelArray::DimensionSpec> ModelArray::definedDimensions = {
    { ModelArray::Dimension::X, { "xdim", 0 } },
    { ModelArray::Dimension::Y, { "ydim", 0 } },
    { ModelArray::Dimension::Z, { "zdim", 1 } },
    { ModelArray::Dimension::XVERTEX, { "xvertex", 1 } }, // defined as x + 1
    { ModelArray::Dimension::YVERTEX, { "yvertex", 1 } }, // defined as y + 1
    { ModelArray::Dimension::XCG, { "x_cg", 1 } },
    { ModelArray::Dimension::YCG, { "y_cg", 1 } },
    // The DG components are also included here to store the names
    { ModelArray::Dimension::DG, { "dg_comp", 1 } },
    { ModelArray::Dimension::DGSTRESS, { "dgstress_comp", 1 } },
    { ModelArray::Dimension::NCOORDS, { "ncoords", 2 } }, // It's a two dimensional model
// clang-format on
};

ModelArray::TypeDimensions ModelArray::typeDimensions = {
    { ModelArray::Type::H,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::VERTEX,
        {
            ModelArray::Dimension::XVERTEX,
            ModelArray::Dimension::YVERTEX,
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
    { ModelArray::Type::DG,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::DGSTRESS,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::CG,
        {
            ModelArray::Dimension::XCG,
            ModelArray::Dimension::YCG,
        } },
};

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::H, "HField" },
    { ModelArray::Type::VERTEX, "VertexField" },
    { ModelArray::Type::U, "UField" },
    { ModelArray::Type::V, "VField" },
    { ModelArray::Type::Z, "ZField" },
    { ModelArray::Type::DG, "DGField" },
    { ModelArray::Type::DGSTRESS, "DGStressField" },
    { ModelArray::Type::CG, "CGField" },
};

ModelArray::ModelArray()
    : ModelArray(Type::H)
{
}

bool ModelArray::hasDoF(const Type type)
{
    return type == Type::DG || type == Type::DGSTRESS || type == Type::VERTEX;
}

ModelArray::SizeMap::SizeMap()
    : m_sizes({ { Type::H, 0 }, { Type::VERTEX, 1 }, { Type::U, 0 }, { Type::V, 0 }, { Type::Z, 0 },
        { Type::DG, 0 }, { Type::DGSTRESS, 0 }, { Type::CG, 1 } })
{
}

ModelArray::DimensionMap::DimensionMap()
    : m_dimensions({
        { Type::H, { 0, 0 } },
        { Type::VERTEX, { 1, 1 } },
        { Type::U, { 0, 0 } },
        { Type::V, { 0, 0 } },
        { Type::Z, { 0, 0, 1 } },
        { Type::DG, { 0, 0 } },
        { Type::DGSTRESS, { 0, 0 } },
        { Type::CG, { 1, 1 } },
    })
{
}
const size_t ModelArray::nCoords = 2;

const std::map<ModelArray::Type, ModelArray::Dimension> ModelArray::componentMap = {
    { Type::DG, Dimension::DG },
    { Type::DGSTRESS, Dimension::DGSTRESS },
    { Type::VERTEX, Dimension::NCOORDS },
};

}
