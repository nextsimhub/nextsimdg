/*!
 * @file ModelArrayRef3.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF3_HPP
#define MODELARRAYREF3_HPP

#include "include/MARBackingStore.hpp"
#include "include/ModelArray.hpp"
#include <map>

namespace Nextsim {

const bool RW = true;
const bool RO = false;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;

template <bool isReadWrite = RO> class ModelArrayRef {
public:
    ModelArrayRef(const std::string& field, MARBackingStore& backingStore)
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
    MARBackingStore store;
    friend MARBackingStore;
};

template <> class ModelArrayRef<RW> {
public:
    ModelArrayRef(const std::string& field, MARBackingStore& backingStore)
    {
        ref = nullptr;
        backingStore.getFieldAddr(field, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    double& operator[](size_t index) const { return ref->operator[](index); }

private:
    ModelArrayReference ref;
    MARBackingStore store;
    friend MARBackingStore;
};
}

#endif /* MODELARRAYREF3_HPP */
