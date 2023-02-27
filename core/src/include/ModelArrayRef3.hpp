/*!
 * @file ModelArrayRef3.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF3_HPP
#define MODELARRAYREF3_HPP

#include "include/MARStore.hpp"
#include "include/ModelArray.hpp"
#include "include/TextTag.hpp"

#include <map>

namespace Nextsim {

const bool RW = true;
const bool RO = false;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;


template <const TextTag& fieldNameTp, bool isReadWrite = RO> class ModelArrayRef {
public:
    ModelArrayRef(MARStore& backingStore)
        : fieldName(fieldNameTp)
        , store(backingStore)
    {
        ref = nullptr;
        store.getFieldAddr(fieldName.text, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    ModelArrayRef(const ModelArrayRef&) = delete;
    ModelArrayRef& operator=(const ModelArrayRef&) = delete;
    const double& operator[](size_t index) const { return ref->operator[](index); }

private:
    ModelArrayConstReference ref;
    MARStore& store;
    TextTag fieldName;
    friend MARStore;
};

template <const TextTag &fieldNameTp> class ModelArrayRef<fieldNameTp, RW> {
public:
    ModelArrayRef(MARStore& backingStore)
        : fieldName(fieldNameTp)
        , store(backingStore)
    {
        ref = nullptr;
        store.getFieldAddr(fieldName.text, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    double& operator[](size_t index) const { return ref->operator[](index); }

private:
    ModelArrayReference ref;
    MARStore& store;
    TextTag fieldName;
    friend MARStore;
};
}

#endif /* MODELARRAYREF3_HPP */
