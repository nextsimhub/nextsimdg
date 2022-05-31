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

// Simple 1d indexing
inline size_t indexr(const size_t* dims, const size_t nDims, size_t i) { return i; }

// Special case for 2d indexing
inline size_t indexr(const size_t* dims, const size_t nDims, size_t i, size_t j)
{
    return i * dims[1] + j;
}

// General case for up to 8d indexing
inline size_t indexr(const size_t* dims, const size_t nDims, size_t i, size_t j, size_t k,
    size_t l = 0, size_t m = 0, size_t n = 0, size_t p = 0, size_t q = 0)
{
    size_t stepLength = 1;
    size_t out = 0;
    switch (nDims) {
    case 8:
        out += stepLength * q;
        stepLength *= dims[7];
    case 7:
        out += stepLength * p;
        stepLength *= dims[6];
    case 6:
        out += stepLength * n;
        stepLength *= dims[5];
    case 5:
        out += stepLength * m;
        stepLength *= dims[4];
    case 4:
        out += stepLength * l;
        stepLength *= dims[3];
    case 3:
        out += stepLength * k;
        stepLength *= dims[2];
    case 2:
        out += stepLength * j;
        stepLength *= dims[1];
    case 1:
    default:
        out += stepLength * i;
        return out;
    }
    return 0;
}

size_t indexr(const ModelArray::Dimensions& loc, const size_t* dims)
{
    size_t loc8[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    for (char i = 0; i < loc.size(); ++i) {
        loc8[i] = loc[i];
    }

    return indexr(
        dims, loc.size(), loc8[0], loc8[1], loc8[2], loc8[3], loc8[4], loc8[5], loc8[6], loc8[7]);
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
    return (*this)(indexr(dimensions().data(), 2, i, j));
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k) const
{
    return (*this)(indexr(dimensions().data(), 3, i, j, k));
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l) const
{
    return (*this)(indexr(dimensions().data(), 4, i, j, k, l));
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
{
    return (*this)(indexr(dimensions().data(), 5, i, j, k, l, m));
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    return (*this)(indexr(dimensions().data(), 6, i, j, k, l, m, n));
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
{
    return (*this)(indexr(dimensions().data(), 7, i, j, k, l, m, n, p));
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
{
    return (*this)(indexr(dimensions().data(), 8, i, j, k, l, m, n, p, q));
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
