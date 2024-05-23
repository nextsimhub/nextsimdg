/*!
 * @file ModelData.cpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

#include <algorithm>
#include <cstdarg>
#include <iterator>
#include <set>
#include <string>
#include <utility>

namespace Nextsim {

ModelArray::SizeMap ModelArray::m_sz;
ModelArray::DimensionMap ModelArray::m_dims;
bool ModelArray::areMapsInvalid = true;

ModelArray::ModelArray(const Type type)
    : type(type)
{
    m_data.resize(std::max(std::size_t { 0 }, m_sz.at(type)), nComponents());
    validateMaps();
}

ModelArray::ModelArray(const ModelArray& orig)
    : ModelArray(orig.type)
{
    setData(orig.m_data);
}

ModelArray& ModelArray::operator=(const ModelArray& orig)
{
    type = orig.type;
    setData(orig.m_data);

    return *this;
}

ModelArray& ModelArray::operator=(const double& fill)
{
    setData(fill);

    return *this;
}

ModelArray ModelArray::operator+(const ModelArray& addend) const
{
    ModelArray result = *this;
    return result += addend;
}

ModelArray ModelArray::operator-(const ModelArray& subtrahend) const
{
    ModelArray result = *this;
    return result -= subtrahend;
}

ModelArray ModelArray::operator*(const ModelArray& multiplier) const
{
    ModelArray result = *this;
    return result *= multiplier;
}

ModelArray ModelArray::operator/(const ModelArray& divisor) const
{
    ModelArray result = *this;
    return result /= divisor;
}

ModelArray ModelArray::operator-() const
{
    ModelArray copy(type);
    copy.m_data = -m_data;
    return copy;
}

ModelArray ModelArray::operator+(const double& x) const
{
    ModelArray result = *this;
    return result += x;
}

ModelArray ModelArray::operator-(const double& x) const
{
    ModelArray result = *this;
    return result -= x;
}

ModelArray ModelArray::operator*(const double& x) const
{
    ModelArray result = *this;
    return result *= x;
}

ModelArray ModelArray::operator/(const double& x) const
{
    ModelArray result = *this;
    return result /= x;
}

ModelArray operator+(const double& x, const ModelArray& y) { return y + x; }

ModelArray operator-(const double& x, const ModelArray& y) { return -(y - x); }

ModelArray operator*(const double& x, const ModelArray& y) { return y * x; }

ModelArray operator/(const double& x, const ModelArray& y)
{
    ModelArray xArray(y.getType());
    xArray.setData(x);
    return xArray /= y;
}

ModelArray ModelArray::max(double max) const
{
    ModelArray maxed = ModelArray(type);
    maxed.m_data.array() = m_data.array().max(max);
    return maxed;
}

ModelArray ModelArray::min(double min) const
{
    ModelArray mined = ModelArray(type);
    mined.m_data.array() = m_data.array().min(min);
    return mined;
}

ModelArray ModelArray::max(const ModelArray& maxArr) const
{
    ModelArray maxed = ModelArray(type);
    maxed.m_data.array() = m_data.array().max(maxArr.m_data);
    return maxed;
}

ModelArray ModelArray::min(const ModelArray& minArr) const
{
    ModelArray mined = ModelArray(type);
    mined.m_data.array() = m_data.array().min(minArr.m_data);
    return mined;
}

ModelArray& ModelArray::clampAbove(double max)
{
    m_data = this->max(max).m_data;
    return *this;
}

ModelArray& ModelArray::clampBelow(double min)
{
    m_data = this->min(min).m_data;
    return *this;
}

ModelArray& ModelArray::clampAbove(const ModelArray& maxArr)
{
    m_data = this->max(maxArr).m_data;
    return *this;
}

ModelArray& ModelArray::clampBelow(const ModelArray& minArr)
{
    m_data = this->min(minArr).m_data;
    return *this;
}

void ModelArray::setData(double value)
{
    resize();
    m_data = value;
}

void ModelArray::setData(const double* pData)
{
    resize();
    auto out = std::copy(pData, pData + m_sz.at(type) * nComponents(), m_data.data());
}

void ModelArray::setData(const DataType& from) { m_data = from; } // setData(from.data()); }

void ModelArray::setData(const ModelArray& from) { setData(from.m_data.data()); }

void ModelArray::setDimensions(Type type, const MultiDim& newDims)
{
    std::vector<Dimension>& dimSpecs = typeDimensions.at(type);
    for (size_t i = 0; i < dimSpecs.size(); ++i) {
        definedDimensions.at(dimSpecs[i]).local_length = newDims[i];
    }
    validateMaps();
}

void ModelArray::setNComponents(std::map<Type, size_t> cMap)
{
    for (auto entry : cMap) {
        setNComponents(entry.first, entry.second);
    }
}

#ifdef USE_MPI
void ModelArray::setDimension(Dimension dim, size_t global_length, size_t local_length, size_t start)
#else
void ModelArray::setDimension(Dimension dim, size_t global_length)
#endif
{
#ifdef USE_MPI
    definedDimensions.at(dim).global_length = global_length;
    definedDimensions.at(dim).local_length = local_length;
    definedDimensions.at(dim).start = start;
#else
    // if MPI is not used then set the local_length to be the same as the global
    definedDimensions.at(dim).global_length = global_length;
    definedDimensions.at(dim).local_length = global_length;
    definedDimensions.at(dim).start = 0;
#endif
    validateMaps();
}

const double& ModelArray::operator[](const MultiDim& loc) const
{
    return (*this)[indexr(this->dimensions().data(), loc)];
}

double& ModelArray::operator[](const MultiDim& dims)
{
    return const_cast<double&>(std::as_const(*this)[dims]);
}

ModelArray::Component ModelArray::components(const MultiDim& loc)
{
    return components(indexr(dimensions().data(), loc));
}

const ModelArray::ConstComponent ModelArray::components(const MultiDim& loc) const
{
    return components(indexr(dimensions().data(), loc));
}

/*!
 * @brief Returns the index for a given set of multi-dimensional location for the given Type.
 *
 * @param type The type to act on.
 * @param loc The multi-dimensional location to return the index for.
 */
size_t ModelArray::indexFromLocation(Type type, const MultiDim& loc)
{
    return indexr(m_dims.at(type).data(), loc);
}

/*!
 * @brief Returns the multi-dimensional location for a given index for the given Type.
 *
 * @param type The type to act on.
 * @param index The index to return the multi-dimensional location for.
 */
ModelArray::MultiDim ModelArray::locationFromIndex(Type type, size_t index)
{
    MultiDim& dims = m_dims.at(type);
    MultiDim loc(dims.size()); // Size to the known number of dimensions
    for (size_t i = 0; i < loc.size(); ++i) {
        size_t theDim = dims[i];
        size_t pos = index % theDim;
        loc[i] = pos;
        index /= theDim;
    }
    return loc;
}

void ModelArray::validateMaps()
{
    m_dims.validate();
    m_sz.validate();
    areMapsInvalid = false;
}

void ModelArray::DimensionMap::validate()
{
    for (auto entry : typeDimensions) {
        Type type = entry.first;
        std::vector<size_t>& dims = m_dimensions[type];
        std::vector<Dimension>& typeDims = entry.second;
        dims.resize(typeDims.size());
        for (size_t i = 0; i < typeDims.size(); ++i) {
            dims[i] = definedDimensions.at(typeDims[i]).local_length;
        }
    }
}

void ModelArray::SizeMap::validate()
{
    for (auto entry : typeDimensions) {
        size_t size = 1;
        std::vector<Dimension>& typeDims = entry.second;
        for (size_t i = 0; i < typeDims.size(); ++i) {
            size *= definedDimensions.at(typeDims[i]).local_length;
        }
        m_sizes.at(entry.first) = size;
    }
}

} /* namespace Nextsim */
