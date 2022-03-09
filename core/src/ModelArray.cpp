/*!
 * @file ModelData.cpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

#include <algorithm>
#include <cstdarg>
#include <set>
#include <utility>

namespace Nextsim {

ModelArray::SizeMap ModelArray::m_sz;
ModelArray::DimensionMap ModelArray::m_dims;

const std::map<ModelArray::Type, std::string> ModelArray::typeNames = {
        {ModelArray::Type::H, "HField"},
        {ModelArray::Type::U, "UField"},
        {ModelArray::Type::V, "VField"},
        {ModelArray::Type::Z, "ZField"},
};

ModelArray::ModelArray(const Type type, const std::string& name)
    : type(type)
    , m_name(name)
{
    if (m_sz.at(type) > 0) {
        m_data.reserve(m_sz.at(type));
    }
}

ModelArray::ModelArray()
    : ModelArray(Type::H, "")
{
}

ModelArray::ModelArray(const ModelArray& orig)
    : ModelArray(orig.type, orig.m_name)
{
    setData(orig.m_data.data());
}

ModelArray& ModelArray::operator=(const ModelArray& orig)
{
    type = orig.type;
    m_name = orig.m_name;
    setData(orig.m_data.data());

    return *this;
}

void ModelArray::setData(const double* pData) { m_data.assign(pData, pData + m_sz.at(type)); }

void ModelArray::setData(const std::vector<double>& from) { setData(from.data()); }

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

const double& ModelArray::operator()(size_t i, size_t j) const
{
    return (*this)(i * m_dims.at(type)[1] + j);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k) const
{
    size_t dim2 = m_dims.at(type)[2];
    return (*this)(i * dim2, j * dim2 + k);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l) const
{
    size_t dim3 = m_dims.at(type)[3];
    return (*this)(i * dim3, j * dim3, k * dim3 + l);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
{
    size_t dim4 = m_dims.at(type)[4];
    return (*this)(i * dim4, j * dim4, k * dim4, l * dim4 + m);
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    size_t dim5 = m_dims.at(type)[5];
    return (*this)(i * dim5, j * dim5, k * dim5, l * dim5, m * dim5 + n);
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
{
    size_t dim6 = m_dims.at(type)[6];
    return (*this)(i * dim6, j * dim6, k * dim6, l * dim6, m * dim6, n * dim6 + p);
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
{
    size_t dim7 = m_dims.at(type)[7];
    return (*this)(i * dim7, j * dim7, k * dim7, l * dim7, m * dim7, n * dim7, p * dim7 + q);
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

} /* namespace Nextsim */
