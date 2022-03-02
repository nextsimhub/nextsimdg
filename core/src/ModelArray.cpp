/*!
 * @file ModelData.cpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"

#include <algorithm>
#include <cstdarg>
#include <utility>

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

const double& ModelArray::operator()(size_t i, size_t j) const { return (*this)(i * m_dims[1] + j); }

const double& ModelArray::operator()(size_t i, size_t j, size_t k) const
{
    return (*this)(i * m_dims[2], j * m_dims[2] + k);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l) const
{
    return (*this)(i * m_dims[3], j * m_dims[3], k * m_dims[3] + l);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
{
    return (*this)(i * m_dims[4], j * m_dims[4], k * m_dims[4], l * m_dims[4] + m);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    return (*this)(i * m_dims[5], j * m_dims[5], k * m_dims[5], l * m_dims[5], m * m_dims[5] + n);
}

const double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
{
    return (*this)(i * m_dims[6], j * m_dims[6], k * m_dims[6], l * m_dims[6], m * m_dims[6],
        n * m_dims[6] + p);
}

const double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
{
    return (*this)(i * m_dims[7], j * m_dims[7], k * m_dims[7], l * m_dims[7], m * m_dims[7],
        n * m_dims[7], p * m_dims[7] + q);
}

double& ModelArray::operator[](const Dimensions& dims) {
    return const_cast<double&>(std::as_const(*this)[dims]);
}
double& ModelArray::operator()(size_t i, size_t j) {
    return const_cast<double&>(std::as_const(*this)(i, j));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k) {
    return const_cast<double&>(std::as_const(*this)(i, j, k));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l) {
    return const_cast<double&>(std::as_const(*this)(i, j, k, l));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) {
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) {
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m, n));
}
double& ModelArray::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) {
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m, n, p));
}
double& ModelArray::operator()(
    size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) {
    return const_cast<double&>(std::as_const(*this)(i, j, k, l, m, n, p, q));
}

} /* namespace Nextsim */
