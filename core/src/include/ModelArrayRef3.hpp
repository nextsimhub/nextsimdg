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

#include <iostream> // FIXME remove me

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
        std::cerr << "RO ModelArrayRef(): field = " << field << ", ref = " << ref << std::endl;
        backingStore.getFieldAddr(field, ref);
        std::cerr << "RO ModelArrayRef(): -> ref = " << ref << std::endl;
    }
    const double& operator[](size_t index) const
    {
        std::cerr << "RO ModelArrayRef[]: ref = " << ref << " ref->size() = " << ref->size() << std::endl;
        return ref->operator[](index);
    }
private:
    ModelArrayConstReference ref;
    friend S;
};

template <typename S> class ModelArrayRef<S, RW> {
public:
    ModelArrayRef(const std::string& field, S& backingStore)
    {
        ref = nullptr;
        std::cerr << "RW ModelArrayRef(): field = " << field << ", ref = " << ref << std::endl;
        backingStore.getFieldAddr(field, ref);
        std::cerr << "RW ModelArrayRef(): -> ref = " << ref << std::endl;
    }
    double& operator[](size_t index) const
    {
        std::cerr << "RW ModelArrayRef[]: ref = " << ref << " ref->size() = " << ref->size() << std::endl;
        return ref->operator[](index);
    }
private:
    ModelArrayReference ref;
    friend S;
};
}

#endif /* MODELARRAYREF3_HPP */
