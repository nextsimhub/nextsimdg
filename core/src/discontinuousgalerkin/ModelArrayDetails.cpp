/*!
 * @file ModelArrayDetails.cpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

namespace Nextsim {
std::map<ModelArray::Dimension, ModelArray::DimensionSpec> ModelArray::definedDimensions = {
    { ModelArray::Dimension::X, { "x", 0 } },
    { ModelArray::Dimension::Y, { "y", 0 } },
    { ModelArray::Dimension::Z, { "z", 1 } },
    { ModelArray::Dimension::XCG, { "x_cg", 1 } },
    { ModelArray::Dimension::YCG, { "y_cg", 1 } },
    // The DG components are also included here to store the names
    { ModelArray::Dimension::DG, { "dg_comp", 1 } },
    { ModelArray::Dimension::DGSTRESS, { "dgstress_comp", 1 } },

};

ModelArray::TypeDimensions ModelArray::typeDimensions = {
    { ModelArray::Type::H,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
        } },
    { ModelArray::Type::U,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
        } },
    { ModelArray::Type::V,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
        } },
    { ModelArray::Type::Z,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Z),
        } },
    { ModelArray::Type::DG,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
        } },
    { ModelArray::Type::DGSTRESS,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::X),
            ModelArray::definedDimensions.at(ModelArray::Dimension::Y),
        } },
    { ModelArray::Type::CG,
        {
            ModelArray::definedDimensions.at(ModelArray::Dimension::XCG),
            ModelArray::definedDimensions.at(ModelArray::Dimension::YCG),
        } },
};

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::H, "HField" },
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

bool ModelArray::hasDoF(const Type type) { return type == Type::DG || type == Type::DGSTRESS; }

ModelArray::SizeMap::SizeMap()
    : m_sizes({ { Type::H, 0 }, { Type::U, 0 }, { Type::V, 0 }, { Type::Z, 0 }, { Type::DG, 0 },
        { Type::DGSTRESS, 0 }, { Type::CG, 1 } })
{
}

ModelArray::DimensionMap::DimensionMap()
    : m_dimensions({
        { Type::H, { 0, 0 } },
        { Type::U, { 0, 0 } },
        { Type::V, { 0, 0 } },
        { Type::Z, { 0, 0, 1 } },
        { Type::DG, { 0, 0 } },
        { Type::DGSTRESS, { 0, 0 } },
        { Type::CG, { 1, 1 } },
    })
{
}

}
