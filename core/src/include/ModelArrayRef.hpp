/*!
 * @file ModelArrayRef.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF_HPP
#define MODELARRAYREF_HPP

#include "include/ModelArray.hpp"
#include "include/TextTag.hpp"

#include <map>
#include "include/ModelArrayReferenceStore.hpp"

namespace Nextsim {

const bool RW = true;
const bool RO = false;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;

/*!
 * @brief A class which provides indirect access to ModelArray.
 *
 * @details Provides access to data from other parts of the model using a
 * TextTag-wrapped string key. The class provides indexing of and access to the
 * referenced ModelArray. Here the returned data are by const references, used
 * for accessing data in a read-only fashion.
 *
 * @tparam fieldName The TextTag containing the name of the field to be referenced.
 * @tparam isReadWrite A boolean which is here false, indicating access to the array is read-only.
 */
template <const TextTag& fieldName, bool isReadWrite = RO> class ModelArrayRef {
public:
    ModelArrayRef(ModelArrayReferenceStore& backingStore)
        : store(backingStore)
    {
        dataReference = nullptr;
        store.getFieldAddr(fieldName.text, dataReference);
    }
    ~ModelArrayRef() { store.removeReference(fieldName.text, dataReference); }
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
    const double& operator[](const ModelArray::MultiDim& dims) { return dataReference->operator[](dims); }
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned.
     *
     * @param index The one dimensional index of the target point.
     */
    const double& operator[](size_t index) const { return dataReference->operator[](index); }
    //! Returns the specified point from a 1 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i) const { return dataReference->operator()(i); }
    //! Returns the specified point from a 2 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j) const { return dataReference->operator()(i, j); }
    //! Returns the specified point from a 3 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j, size_t k) const
    {
        return dataReference->operator()(i, j, k);
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
        return dataReference->zIndexAndLayer(hIndex, layer);
    }

    //! Direct access top the underlying data array.
    const ModelArray& data() const { return *dataReference; }
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

    //! Returns a ModelArray containing the per-element sum of the
    //! object and the provided ModelArray.
    ModelArray operator+(double addend) const { return data() + addend; }
    //! Returns a ModelArray containing the per-element difference between the
    //! object and the provided ModelArray.
    ModelArray operator-(double subtrahend) const { return data() - subtrahend; }
    //! Returns a ModelArray containing the per-element product of the
    //! object and the provided ModelArray.
    ModelArray operator*(double multiplier) const { return data() * multiplier; }
    //! Returns a ModelArray containing the per-element ratio between the
    //! object and the provided ModelArray.
    ModelArray operator/(double divisor) const { return data() / divisor; }

private:
    ModelArrayConstReference dataReference;
    ModelArrayReferenceStore& store;
    friend ModelArrayReferenceStore;
};

/*!
 * @brief A class which provides indirect access to ModelArray.
 *
 * @details Provides access to data from other parts of the model using a
 * TextTag-wrapped string key. The class provides indexing of and access to the
 * referenced ModelArray. Here the returned data are by non-const references, used
 * for accessing data in a read-write fashion.
 *
 * @tparam fieldName The TextTag containing the name of the field to be referenced.
 * @tparam isReadWrite A boolean which is here true, indicating access to the array is read-write.
 */
template <const TextTag& fieldName> class ModelArrayRef<fieldName, RW> {
public:
    ModelArrayRef(ModelArrayReferenceStore& backingStore)
        : store(backingStore)
    {
        dataReference = nullptr;
        store.getFieldAddr(fieldName.text, dataReference);
    }
    ~ModelArrayRef() { store.removeReference(fieldName.text, dataReference); }
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
    double& operator[](const ModelArray::MultiDim& dims) { return dataReference->operator[](dims); }
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned.
     *
     * @param index The one dimensional index of the target point.
     */
    double& operator[](size_t index) const { return dataReference->operator[](index); }
    //! Returns the specified point from a 1 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i) const { return dataReference->operator()(i); }
    //! Returns the specified point from a 2 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j) const { return dataReference->operator()(i, j); }
    //! Returns the specified point from a 3 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k) const { return dataReference->operator()(i, j, k); }

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
        return dataReference->zIndexAndLayer(hIndex, layer);
    }

    //! Direct access top the underlying data array.
    ModelArray& data() const { return *dataReference; }
    //! Cast the reference class to a real reference to the referenced ModelArray.
    operator ModelArray&() const { return data(); }

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

    //! Returns a ModelArray containing the per-element sum of the
    //! object and the provided ModelArray.
    ModelArray operator+(double addend) const { return data() + addend; }
    //! Returns a ModelArray containing the per-element difference between the
    //! object and the provided ModelArray.
    ModelArray operator-(double subtrahend) const { return data() - subtrahend; }
    //! Returns a ModelArray containing the per-element product of the
    //! object and the provided ModelArray.
    ModelArray operator*(double multiplier) const { return data() * multiplier; }
    //! Returns a ModelArray containing the per-element ratio between the
    //! object and the provided ModelArray.
    ModelArray operator/(double divisor) const { return data() / divisor; }

private:
    ModelArrayReference dataReference;
    ModelArrayReferenceStore& store;
    friend ModelArrayReferenceStore;
};
}

#endif /* MODELARRAYREF_HPP */
