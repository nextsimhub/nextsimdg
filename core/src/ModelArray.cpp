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

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
    { ModelArray::Type::H, "HField" },
    { ModelArray::Type::U, "UField" },
    { ModelArray::Type::V, "VField" },
    { ModelArray::Type::Z, "ZField" },
    { ModelArray::Type::DG, "DGHField--DO-NOT-USE--" },
};

ModelArray::ModelArray(const Type type, const std::string& name)
    : type(type)
    , m_name(name)
{
    m_data.resize(std::max(std::size_t { 0 }, m_sz.at(type)), nComponents());
}

ModelArray::ModelArray()
    : ModelArray(Type::H, "")
{
}

ModelArray::ModelArray(const ModelArray& orig)
    : ModelArray(orig.type, orig.m_name)
{
    setData(orig.m_data);
}

ModelArray& ModelArray::operator=(const ModelArray& orig)
{
    type = orig.type;
    m_name = orig.m_name;
    setData(orig.m_data);

    return *this;
}

ModelArray& ModelArray::operator=(const double& fill)
{
    m_name = std::to_string(fill);
    setData(fill);

    return *this;
}

ModelArray ModelArray::operator+(const ModelArray& addend) const
{
    ModelArray result(type, m_name + "+" + addend.m_name);

    for (size_t i = 0; i < m_sz.at(type); ++i) {
        result[i] = (*this)[i] + addend[i];
    }

    return result;
}

ModelArray ModelArray::operator-(const ModelArray& subtrahend) const
{
    ModelArray result(type, m_name + "-" + subtrahend.m_name);

    for (size_t i = 0; i < m_sz.at(type); ++i) {
        result[i] = (*this)[i] - subtrahend[i];
    }

    return result;
}

ModelArray ModelArray::operator*(const ModelArray& multiplier) const
{
    ModelArray result(type, m_name + "*" + multiplier.m_name);

    for (size_t i = 0; i < m_sz.at(type); ++i) {
        result[i] = (*this)[i] * multiplier[i];
    }

    return result;
}

ModelArray ModelArray::operator/(const ModelArray& divisor) const
{
    ModelArray result(type, m_name + "/" + divisor.m_name);

    for (size_t i = 0; i < m_sz.at(type); ++i) {
        result[i] = (*this)[i] / divisor[i];
    }

    return result;
}

ModelArray ModelArray::operator-() const
{
    ModelArray copy(type, std::string("-") + m_name);
    copy.m_data = -m_data;
    return copy;
}

ModelArray ModelArray::operator+(const double& x) const
{
    ModelArray copy(type, m_name + "+" + std::to_string(x));
    copy.m_data = m_data + x;
    return copy;
}

ModelArray ModelArray::operator-(const double& x) const
{
    ModelArray copy(type, m_name + "-" + std::to_string(x));
    copy.m_data = m_data - x;
    return copy;
}

ModelArray ModelArray::operator*(const double& x) const
{
    ModelArray copy(type, m_name + "*" + std::to_string(x));
    copy.m_data = m_data * x;
    return copy;
}

ModelArray ModelArray::operator/(const double& x) const
{
    ModelArray copy(type, m_name + "/" + std::to_string(x));
    copy.m_data = m_data / x;
    return copy;
}

ModelArray operator+(const double& x, const ModelArray& y) { return y + x; }

ModelArray operator-(const double& x, const ModelArray& y) { return -(y - x); }

ModelArray operator*(const double& x, const ModelArray& y) { return y * x; }

ModelArray operator/(const double& x, const ModelArray& y)
{
    ModelArray xArray(y.getType(), std::to_string(x));
    xArray.setData(x);
    return xArray / y;
}

void ModelArray::setData(double value)
{
    resize();
    m_data = value;
}

void ModelArray::setData(const double* pData)
{
    resize();
    auto out = std::copy(pData, pData + m_sz.at(type), m_data.data());
}

void ModelArray::setData(const DataType& from) { setData(from.data()); }

void ModelArray::setData(const ModelArray& from) { setData(from.m_data.data()); }

void ModelArray::setDimensions(Type type, const Dimensions& newDims)
{
    size_t newSize = 1;
    for (size_t dimLen : newDims) {
        newSize *= dimLen;
    }
    m_dims.at(type).clear();
    std::copy(newDims.begin(), newDims.end(), std::back_inserter(m_dims.at(type)));
    m_sz.at(type) = newSize;
}

const double& ModelArray::operator[](const Dimensions& dims) const
{
    switch (dims.size()) {
    case (1):
        return (*this)(dims[0]);
    case (2):
        return (*this)(dims[0], dims[1]);
    case (3):
        return (*this)(dims[0], dims[1], dims[2]);
    case (4):
        return (*this)(dims[0], dims[1], dims[2], dims[3]);
    case (5):
        return (*this)(dims[0], dims[1], dims[2], dims[3], dims[4]);
    case (6):
        return (*this)(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5]);
    case (7):
        return (*this)(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6]);
    default:
        return (*this)(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6], dims[7]);
    }
}

double& ModelArray::operator[](const Dimensions& dims)
{
    return const_cast<double&>(std::as_const(*this)[dims]);
}

ModelArray::Component ModelArray::components(const Dimensions& loc)
{
    return components(indexr(loc, dimensions().data()));
}
} /* namespace Nextsim */
