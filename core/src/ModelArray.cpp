/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

#include "include/ModelArraySlice.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <utility>

namespace Nextsim {

ModelArray::SizeMap ModelArray::m_sz;
ModelArray::DimensionMap ModelArray::m_dims;
bool ModelArray::areMapsInvalid = true;

ModelArray::ModelArray(const Type type, const std::pair<double, double>& bounds)
    : type(type)
{
    m_data.resize(std::max(std::size_t { 0 }, m_sz.at(type)), nComponents());
    validateMaps();
    setLimits(bounds.first, bounds.second);
}

ModelArray::ModelArray(const ModelArray& orig)
    : ModelArray(orig.type)
{
    setData(orig.m_data);
}

ModelArray& ModelArray::operator=(const ModelArray& orig)
{
    if (orig.nComponents() == 1 && nComponents() != 1) {
        component(0) = orig.data();
    } else {
        type = orig.type;
        setData(orig.m_data);
    }

    return *this;
}

ModelArray& ModelArray::operator=(const double& fill)
{
    setData(fill);

    return *this;
}

ModelArray& ModelArray::operator=(const ModelArraySlice& mas)
{
    return mas.copyToModelArray(*this);
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

ModelArray ModelArray::sqrt()
{
    ModelArray sqrted = ModelArray(type);
    sqrted.m_data.array() = m_data.array().sqrt();
    return sqrted;
}

ModelArray ModelArray::sin()
{
    auto sined = ModelArray(type);
    sined.m_data.array() = m_data.array().sin();
    return sined;
}

ModelArray ModelArray::cos()
{
    auto cosed = ModelArray(type);
    cosed.m_data.array() = m_data.array().cos();
    return cosed;
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
        definedDimensions.at(dimSpecs[i]).localLength = newDims[i];
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
void ModelArray::setDimension(Dimension dim, size_t globalLength, size_t localLength, size_t start)
{
    definedDimensions.at(dim).setLengths(globalLength, localLength, start);
#else
void ModelArray::setDimension(Dimension dim, size_t globalLength)
{
    definedDimensions.at(dim).setLengths(globalLength);
#endif
    validateMaps();
}

const double& ModelArray::operator[](const MultiDim& loc) const
{
    return (*this)[indexr(this->dimensions(), loc)];
}

double& ModelArray::operator[](const MultiDim& dims)
{
    return const_cast<double&>(std::as_const(*this)[dims]);
}

ModelArraySlice ModelArray::operator[](const Slice& slice) { return ModelArraySlice(*this, slice); }

ConstModelArraySlice ModelArray::operator[](const Slice& slice) const
{
    return ConstModelArraySlice(*this, slice);
}

ModelArray::Components ModelArray::components(const MultiDim& loc)
{
    return components(indexr(dimensions(), loc));
}

const ModelArray::ConstComponents ModelArray::components(const MultiDim& loc) const
{
    return components(indexr(dimensions(), loc));
}

/*!
 * @brief Returns the index for a given set of multi-dimensional location for the given Type.
 *
 * @param type The type to act on.
 * @param loc The multi-dimensional location to return the index for.
 */
size_t ModelArray::indexFromLocation(Type type, const MultiDim& loc)
{
    return indexr(m_dims.at(type), loc);
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

void ModelArray::setLimits(const double lower, const double upper)
{
    lowerPhysicalLimit = lower;
    upperPhysicalLimit = upper;
    fillValue = (lowerPhysicalLimit + upperPhysicalLimit) * 0.5;
}

void ModelArray::checkLimits(const ModelArray& mask) const
{
    // Mask the data with the land mask
    const auto masked = (mask.data() == 1).select(m_data.col(0), fillValue);

    // Check first for NaNs. The code is different for the bounds check, because Eigen doesn't
    // return an index for NaN-checking.
    if (masked.isNaN().any())
        throw std::runtime_error("Field contains NaN.");

    /* Now we check the bounds and set the array index (i) and value if we're out of bounds.
     * Here, we need to check if the values are _outside_ the bounds, and if they are, then we ask
     * Eigen to find the offending value and its location. We then proceed to throw an error.
     * This also means that using '<' and '>' in the checks here is consistent with checking if the
     * value is in min <= value <= max.
     */
    size_t i;
    double value;
    if (masked.minCoeff() < lowerPhysicalLimit) {
        value = masked.col(0).minCoeff(&i);
    } else if (masked.maxCoeff() > upperPhysicalLimit) {
        value = masked.col(0).maxCoeff(&i);
    } else {
        return;
    }

    /* If we haven't returned (or thrown an exception) by now, we have an error in the field, and
     * Eigen has found that this is at index i.
     */
    const std::vector<size_t> loc = locationFromIndex(type, i);
    std::string locStr = "[";
    for (const size_t& l : loc)
        locStr += std::to_string(l) + ",";
    locStr.pop_back();
    locStr.push_back(']');

    throw std::runtime_error("Field contains out-of-bounds value(s), " + std::to_string(value)
        + " not in [" + std::to_string(lowerPhysicalLimit) + ","
        + std::to_string(upperPhysicalLimit) + "]. Error at index " + locStr + ".");
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
            dims[i] = definedDimensions.at(typeDims[i]).localLength;
        }
    }
}

void ModelArray::SizeMap::validate()
{
    for (auto entry : typeDimensions) {
        size_t size = 1;
        std::vector<Dimension>& typeDims = entry.second;
        for (size_t i = 0; i < typeDims.size(); ++i) {
            size *= definedDimensions.at(typeDims[i]).localLength;
        }
        m_sizes.at(entry.first) = size;
    }
}

} /* namespace Nextsim */
