/*!
 * @file ModelData.cpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelData.hpp"

#include <algorithm>
#include <cstdarg>

namespace Nextsim {

size_t ModelData::m_sz;
ModelData::Dimensions ModelData::m_dims;

ModelData::ModelData(const std::string& name)
    : m_name(name)
{
    if (m_sz > 0) {
        m_data.reserve(m_sz);
    }
}

ModelData::ModelData()
    : ModelData("")
{}

void ModelData::setData(double* pData)
{
    m_data.assign(pData, pData + m_sz);
}

void ModelData::setData(const std::vector<double>& from) {
    m_data.clear();
    std::copy(from.begin(), from.end(), std::back_inserter(m_data));
}

void ModelData::setDimensions(const Dimensions& newDims)
{
    size_t newSize = 1;
    for (size_t dimLen: newDims) {
        newSize *= dimLen;
    }
    m_dims.clear();
    std::copy(newDims.begin(), newDims.end(), std::back_inserter(m_dims));
    m_sz = newSize;
}

double& ModelData::operator[](const Dimensions& dims)
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

double& ModelData::operator()(size_t i, size_t j) { return (*this)(i * m_dims[1] + j); }

double& ModelData::operator()(size_t i, size_t j, size_t k)
{
    return (*this)(i * m_dims[2], j * m_dims[2] + k);
}

double& ModelData::operator()(size_t i, size_t j, size_t k, size_t l)
{
    return (*this)(i * m_dims[3], j * m_dims[3], k * m_dims[3] + l);
}

double& ModelData::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
    return (*this)(i * m_dims[4], j * m_dims[4], k * m_dims[4], l * m_dims[4] + m);
}

double& ModelData::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
    return (*this)(i * m_dims[5], j * m_dims[5], k * m_dims[5], l * m_dims[5], m * m_dims[5] + n);
}

double& ModelData::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p)
{
    return (*this)(i * m_dims[6], j * m_dims[6], k * m_dims[6], l * m_dims[6], m * m_dims[6],
        n * m_dims[6] + p);
}

double& ModelData::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q)
{
    return (*this)(i * m_dims[7], j * m_dims[7], k * m_dims[7], l * m_dims[7], m * m_dims[7],
        n * m_dims[7], p * m_dims[7] + q);
}

} /* namespace Nextsim */
