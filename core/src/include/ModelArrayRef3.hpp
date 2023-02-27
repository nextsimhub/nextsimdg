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

struct TextTag {
    operator std::string() const { return std::string(text); };
    const char* text;
};

template <const TextTag &fieldNameTp, bool isReadWrite = RO> class ModelArrayRef {
public:
    ModelArrayRef(MARBackingStore& backingStore)
        : fieldName(fieldNameTp)
    {
        ref = nullptr;
        backingStore.getFieldAddr(fieldName.text, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    ModelArrayRef(const ModelArrayRef&) = delete;
    ModelArrayRef& operator=(const ModelArrayRef&) = delete;
    const double& operator[](size_t index) const { return ref->operator[](index); }

private:
    ModelArrayConstReference ref;
    MARBackingStore store;
    TextTag fieldName;
    friend MARBackingStore;
};

template <const TextTag &fieldNameTp> class ModelArrayRef<fieldNameTp, RW> {
public:
    ModelArrayRef(MARBackingStore& backingStore)
        : fieldName(fieldNameTp)
    {
        ref = nullptr;
        backingStore.getFieldAddr(fieldName.text, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    double& operator[](size_t index) const { return ref->operator[](index); }

private:
    TextTag fieldName;
    ModelArrayReference ref;
    MARBackingStore store;
    friend MARBackingStore;
};
}

#endif /* MODELARRAYREF3_HPP */
