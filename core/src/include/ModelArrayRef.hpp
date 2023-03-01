/*!
 * @file ModelArrayRef.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF_HPP
#define MODELARRAYREF_HPP

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
    /*!
     * @brief Returns the data at the indices.
     *
     * @details The argument is a list of dimension indices (actually a
     * std::vector<size_t>). The number of dimensions provided can be lower
     * than that of the ModelArray type. If the object holds discontinuous
     * Galerkin components, only the cell averaged value is returned.
     *
     * @param dims The indices of the target point.
     */
    const double& operator[](const ModelArray::MultiDim& dims) { return ref->operator[](dims); }
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned.
     *
     * @param index The one dimensional index of the target point.
     */
    const double& operator[](size_t index) const { return ref->operator[](index); }
    //! Returns the specified point from a 1 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i) const { return ref->operator()(i); }
    //! Returns the specified point from a 2 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j) const { return ref->operator()(i, j); }
    //! Returns the specified point from a 3 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j, size_t k) const
    {
        return ref->operator()(i, j, k);
    }

    /*!
     * @brief Special access function for ZFields.
     *
     * @detail Index the referenced ZField using an index from an HField of the
     * same horizontal extent and a layer index for the final dimension.
     *
     * @param hIndex the equivalent positional index in an HField array
     * @param layer the vertical layer to be accessed
     */
    const double& zIndexAndLayer(size_t hIndex, size_t layer)
    {
        return ref->zIndexAndLayer(hIndex, layer);
    }

    //! Direct access top the underlying data array.
    const ModelArray& data() const { return *ref; }
    //! Cast the reference class to a real reference to the referenced ModelArray.
    operator const ModelArray&() const { return data(); }

    //! Returns a ModelArray containing the per-element sum of the
    //! object and the provided ModelArray.
    ModelArray operator+(const ModelArray& addend) const { return data() + addend; }
    //! Returns a ModelArray containing the per-element difference between the
    //! object and the provided ModelArray.
    ModelArray operator-(const ModelArray& subtrahend) const { return data() - subtrahend; }
    //! Returns a ModelArray containing the per-element product of the
    //! object and the provided ModelArray.
    ModelArray operator*(const ModelArray& multiplier) const { return data() * multiplier; }
    //! Returns a ModelArray containing the per-element ratio between the
    //! object and the provided ModelArray.
    ModelArray operator/(const ModelArray& divisor) const { return data() / divisor; }

private:
    ModelArrayConstReference ref;
    MARStore& store;
    TextTag fieldName;
    friend MARStore;
};

template <const TextTag& fieldNameTp> class ModelArrayRef<fieldNameTp, RW> {
public:
    ModelArrayRef(MARStore& backingStore)
        : fieldName(fieldNameTp)
        , store(backingStore)
    {
        ref = nullptr;
        store.getFieldAddr(fieldName.text, ref);
    }
    ~ModelArrayRef() { store.removeReference(&ref); }
    /*!
     * @brief Returns the data at the indices.
     *
     * @details The argument is a list of dimension indices (actually a
     * std::vector<size_t>). The number of dimensions provided can be lower
     * than that of the ModelArray type. If the object holds discontinuous
     * Galerkin components, only the cell averaged value is returned.
     *
     * @param dims The indices of the target point.
     */
    double& operator[](const ModelArray::MultiDim& dims) { return ref->operator[](dims); }
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned.
     *
     * @param index The one dimensional index of the target point.
     */
    double& operator[](size_t index) const { return ref->operator[](index); }
    //! Returns the specified point from a 1 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i) const { return ref->operator()(i); }
    //! Returns the specified point from a 2 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j) const { return ref->operator()(i, j); }
    //! Returns the specified point from a 3 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k) const { return ref->operator()(i, j, k); }

    /*!
     * @brief Special access function for ZFields.
     *
     * @detail Index the referenced ZField using an index from an HField of the
     * same horizontal extent and a layer index for the final dimension.
     *
     * @param hIndex the equivalent positional index in an HField array
     * @param layer the vertical layer to be accessed
     */
    double& zIndexAndLayer(size_t hIndex, size_t layer)
    {
        return ref->zIndexAndLayer(hIndex, layer);
    }

    //! Direct access top the underlying data array.
    ModelArray& data() const { return *ref; }
    //! Cast the reference class to a real reference to the referenced ModelArray.
    operator ModelArray&() const { return data(); }

private:
    ModelArrayReference ref;
    MARStore& store;
    TextTag fieldName;
    friend MARStore;
};
}

#endif /* MODELARRAYREF_HPP */
