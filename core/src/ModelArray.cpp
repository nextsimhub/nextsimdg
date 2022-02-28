/*!
 * @file ModelData.cpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

#include <algorithm>
#include <cstdarg>

namespace Nextsim {

size_t ModelArray::m_sz;
ModelArray::Dimensions ModelArray::m_dims;

ModelArray::ModelArray(const std::string& name)
    : m_name(name)
{
    if (m_sz > 0) {
        m_data.reserve(m_sz);
    }
}

ModelArray::ModelArray()
    : ModelArray("")
{}

ModelArray::ModelArray(const ModelArray& orig)
    : ModelArray(orig.m_name)
{
    setData(orig.m_data.data());
}

ModelArray& ModelArray::operator=(const ModelArray& orig)
{
    m_name = orig.m_name;
    setData(orig.m_data.data());

    return *this;
}

void ModelArray::setData(const double* pData)
{
    m_data.assign(pData, pData + m_sz);
}

void ModelArray::setData(const std::vector<double>& from) {
    setData(from.data());
}

void ModelArray::setDimensions(const Dimensions& newDims)
{
    size_t newSize = 1;
    for (size_t dimLen: newDims) {
        newSize *= dimLen;
    }
    m_dims.clear();
    std::copy(newDims.begin(), newDims.end(), std::back_inserter(m_dims));
    m_sz = newSize;
}

double& ModelArray::operator[](const Dimensions& dims)
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

double& ModelArray::operator()(size_t i, size_t j) { return (*this)(i * m_dims[1] + j); }

double& ModelArray::operator()(size_t i, size_t j, size_t k)
{
    return (*this)(i * m_dims[2], j * m_dims[2] + k);
}

double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l)
{
    return (*this)(i * m_dims[3], j * m_dims[3], k * m_dims[3] + l);
}

double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
    return (*this)(i * m_dims[4], j * m_dims[4], k * m_dims[4], l * m_dims[4] + m);
}

double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
    return (*this)(i * m_dims[5], j * m_dims[5], k * m_dims[5], l * m_dims[5], m * m_dims[5] + n);
}

double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p)
{
    return (*this)(i * m_dims[6], j * m_dims[6], k * m_dims[6], l * m_dims[6], m * m_dims[6],
        n * m_dims[6] + p);
}

double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q)
{
    return (*this)(i * m_dims[7], j * m_dims[7], k * m_dims[7], l * m_dims[7], m * m_dims[7],
        n * m_dims[7], p * m_dims[7] + q);
}

} /* namespace Nextsim */
