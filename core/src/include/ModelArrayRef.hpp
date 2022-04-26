/*!
 * @file ModelArrayRef.hpp
 *
 * @date Apr 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF_HPP
#define MODELARRAYREF_HPP

#include "ModelArray.hpp"
#include "ModelComponent.hpp"
namespace Nextsim {
const bool RW = true;
const bool RO = false;

template <auto autoType, bool access = RO> class ModelArrayRef {
public:
    const double& operator[](const ModelArray::Dimensions& dims)
    {
        return ModelComponent::getConstArray<autoType>()->operator[](dims);
    }
    const double& operator[](size_t index) const
    {
        return ModelComponent::getConstArray<autoType>()->operator[](index);
    }
    const double& operator()(size_t i) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i);
    }
    const double& operator()(size_t i, size_t j) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j);
    }
    const double& operator()(size_t i, size_t j, size_t k) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j, k);
    }
    const double& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j, k, l);
    }
    const double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j, k, l, m);
    }
    const double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j, k, l, m, n);
    }
    const double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j, k, l, m, n, p);
    }
    const double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
    {
        return ModelComponent::getConstArray<autoType>()->operator()(i, j, k, l, m, n, p, q);
    }

    const double& zIndexAndLayer(size_t hIndex, size_t layer)
    {
        return ModelComponent::getConstArray<autoType>()->zIndexAndLayer(hIndex, layer);
    }

    const ModelArray& data() const { return *ModelComponent::getConstArray<autoType>(); }
    operator const ModelArray&() const { return data(); }
};

template <ModelComponent::SharedArray sh> class ModelArrayRef<sh, RW> {
public:
    double& operator[](const ModelArray::Dimensions& dims)
    {
        return ModelComponent::getArray<sh>()->operator[](dims);
    }
    double& operator[](size_t index) const
    {
        return ModelComponent::getArray<sh>()->operator[](index);
    }
    double& operator()(size_t i) const { return ModelComponent::getArray<sh>()->operator()(i); }
    double& operator()(size_t i, size_t j) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j);
    }
    double& operator()(size_t i, size_t j, size_t k) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j, k);
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j, k, l);
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j, k, l, m);
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j, k, l, m, n);
    }
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j, k, l, m, n, p);
    }
    double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
    {
        return ModelComponent::getArray<sh>()->operator()(i, j, k, l, m, n, p, q);
    }

    double& zIndexAndLayer(size_t hIndex, size_t layer)
    {
        return ModelComponent::getArray<sh>()->zIndexAndLayer(hIndex, layer);
    }

    ModelArray& data() const { return *ModelComponent::getArray<sh>(); }
    operator ModelArray&() const { return data(); }

};
}
#endif /* MODELARRAYREF_HPP */
