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

inline size_t indexr(size_t i, size_t j, const size_t* d) { return i * d[1] + j; }

inline size_t indexr(size_t i, size_t j, size_t k, const size_t* d)
{
    size_t dLast = d[2];
    return indexr(i * dLast, j * dLast + k, d);
}

inline size_t indexr(size_t i, size_t j, size_t k, size_t l, const size_t* d)
{
    size_t dLast = d[3];
    return indexr(i * dLast, j * dLast, k * dLast + l, d);
}

inline size_t indexr(size_t i, size_t j, size_t k, size_t l, size_t m, const size_t* d)
{
    size_t dLast = d[4];
    return indexr(i * dLast, j * dLast, k * dLast, l * dLast + m, d);
}

inline size_t indexr(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, const size_t* d)
{
    size_t dLast = d[5];
    return indexr(i * dLast, j * dLast, k * dLast, l * dLast, m * dLast + n, d);
}

inline size_t indexr(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, const size_t* d)
{
    size_t dLast = d[6];
    return indexr(i * dLast, j * dLast, k * dLast, l * dLast, m * dLast, n * dLast + p, d);
}

inline size_t indexr(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q, const size_t* d)
{
    size_t dLast = d[6];
    return indexr(
        i * dLast, j * dLast, k * dLast, l * dLast, m * dLast, n * dLast, p * dLast + q, d);
}

size_t indexr(const ModelArray::Dimensions& loc, const size_t* dims)
{
    switch (loc.size()) {
    case (1):
        return loc[0];
    case (2):
        return indexr(loc[0], loc[1], dims);
    case (3):
        return indexr(loc[0], loc[1], loc[2], dims);
    case (4):
        return indexr(loc[0], loc[1], loc[2], loc[3], dims);
    case (5):
        return indexr(loc[0], loc[1], loc[2], loc[3], loc[4], dims);
    case (6):
        return indexr(loc[0], loc[1], loc[2], loc[3], loc[4], loc[5], dims);
    case (7):
        return indexr(loc[0], loc[1], loc[2], loc[3], loc[4], loc[5], loc[6], dims);
    default:
        return indexr(loc[0], loc[1], loc[2], loc[3], loc[4], loc[5], loc[6], loc[7], dims);
    }
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

const double& ModelArray::operator()(size_t i, size_t j) const
{
    return (*this)(indexr(i, j, dimensions().data()));
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k) const
{
    return (*this)(indexr(i, j, k, dimensions().data()));
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l) const
{
    return (*this)(indexr(i, j, k, l, dimensions().data()));
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
{
    return (*this)(indexr(i, j, k, l, m, dimensions().data()));
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    return (*this)(indexr(i, j, k, l, m, n, dimensions().data()));
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
{
    return (*this)(indexr(i, j, k, l, m, n, p, dimensions().data()));
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
{
    return (*this)(indexr(i, j, k, l, m, n, p, q, dimensions().data()));
}

double& ModelArray::operator[](const Dimensions& dims)
{
    return const_cast<double&>(std::as_const(*this)[dims]);
}
double& ModelArray::operator()(size_t i, size_t j)
{
    return const_cast<double&>(std::as_const(*this)(i, j));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k)
{
    return const_cast<double&>(std::as_const(*this)(i, j, k));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l)
{
    return const_cast<double&>(std::as_const(*this)(i, j, k, l));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m, n));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p)
{
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m, n, p));
}
double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q)
{
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m, n, p, q));
}

ModelArray::Component ModelArray::components(const Dimensions& loc)
{
    return components(indexr(loc, dimensions().data()));
}
} /* namespace Nextsim */
