/*!
 * @file ModelArrayDetails.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

namespace Nextsim {
std::map<ModelArray::Dimension, ModelArray::DimensionSpec> ModelArray::definedDimensions = {
    { ModelArray::Dimension::X, { "x", "x", 1, 1, 0 } },
    { ModelArray::Dimension::Y, { "y", "y", 1, 1, 0 } },
    { ModelArray::Dimension::Z, { "z", "z", 1, 1, 0 } },
    { ModelArray::Dimension::U, { "u", "u", 1, 1, 0 } },
};

ModelArray::TypeDimensions ModelArray::typeDimensions = {
    { ModelArray::Type::ONED,
        {
            ModelArray::Dimension::X,
        } },
    { ModelArray::Type::TWOD,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
        } },
    { ModelArray::Type::ZUFIELD,
        {
            ModelArray::Dimension::Z,
            ModelArray::Dimension::U,
        } },
    { ModelArray::Type::THREED,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
            ModelArray::Dimension::Z,
        } },
    { ModelArray::Type::FOURD,
        {
            ModelArray::Dimension::X,
            ModelArray::Dimension::Y,
            ModelArray::Dimension::Z,
            ModelArray::Dimension::U,
        } },
};

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::ONED, "OneDField" },
    { ModelArray::Type::TWOD, "TwoDField" },
    { ModelArray::Type::ZUFIELD, "TwoDField" },
    { ModelArray::Type::THREED, "ThreeDField" },
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
          { Type::ZUFIELD, 1 },
          { Type::THREED, 1 },
          { Type::FOURD, 1 },
      })
{
}

ModelArray::DimensionMap::DimensionMap()
    : m_dimensions({
          { Type::ONED, { 1 } },
          { Type::TWOD, { 1, 1 } },
          { Type::ZUFIELD, { 1, 1 } },
          { Type::THREED, { 1, 1, 1 } },
          { Type::FOURD, { 1, 1, 1, 1 } },
      })
{
}

const std::map<ModelArray::Type, ModelArray::Dimension> ModelArray::componentMap = {};

}
