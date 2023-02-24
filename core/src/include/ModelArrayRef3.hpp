/*!
 * @file ModelArrayRef3.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF3_HPP
#define MODELARRAYREF3_HPP

#include "ModelArray.hpp"
#include <map>

namespace Nextsim {

const bool RW = true;
const bool RO = false;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;

template <typename S, bool isReadWrite = RO> class ModelArrayRef {
public:
    ModelArrayRef(const std::string& field, S& backingStore)
    {
        ref = nullptr;
        backingStore.getFieldAddr(field, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    ModelArrayRef(const ModelArrayRef&) = delete;
    ModelArrayRef& operator=(const ModelArrayRef&) = delete;
    const double& operator[](size_t index) const { return ref->operator[](index); }

private:
    ModelArrayConstReference ref;
    S store;
    friend S;
};

template <typename S> class ModelArrayRef<S, RW> {
public:
    ModelArrayRef(const std::string& field, S& backingStore)
    {
        ref = nullptr;
        backingStore.getFieldAddr(field, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    double& operator[](size_t index) const { return ref->operator[](index); }

private:
    ModelArrayReference ref;
    S store;
    friend S;
};
}

#endif /* MODELARRAYREF3_HPP */
