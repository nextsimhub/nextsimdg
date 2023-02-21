/*!
 * @file ModelArrayRef.hpp
 *
 * @date Apr 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYREF2_HPP
#define MODELARRAYREF2_HPP

#include "ModelArray.hpp"
#include <vector>

namespace Nextsim {
const bool RW = true;
const bool RO = false;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;
typedef std::vector<ModelArrayReference> MARBackingStore;
typedef std::vector<ModelArrayConstReference> MARConstBackingStore;

/*!
 * @brief A class which provides indirect access to ModelArray.
 *
 * @details Uses the ModelComponent registry of ModelArrays to access data from
 * other ModelComponents. The class provides indexing, casting and access (by
 * reference) to the pointed to ModelArray. Here the returned data are by const
 * references, used for  ProtectedArray and SharedArray in read-only (RO) mode.
 *
 * @tparam autoType The ModelComponent::ProtectedArray of
 *           ModelComponent::SharedArray to be referenced.
 * @tparam access Whether access to the data is read only (RO = false) or
 *           read-write (RW = true)
 */
template <auto arrayName, typename S, bool access = RO> class ModelArrayRef {
public:
    ModelArrayRef(const S& backingStoreIn)
        : backingStore(backingStoreIn)
    {
    }

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
    const double& operator[](const ModelArray::MultiDim& dims)
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator[](dims);
    }
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned.
     *
     * @param index The one dimensional index of the target point.
     */
    const double& operator[](size_t index) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator[](index);
    }
    //! Returns the specified point from a 1 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i);
    }
    //! Returns the specified point from a 2 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j);
    }
    //! Returns the specified point from a 3 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j, size_t k) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k);
    }
    //! Returns the specified point from a 4 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l);
    }
    //! Returns the specified point from a 5 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m);
    }
    //! Returns the specified point from a 6 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m, n);
    }
    //! Returns the specified point from a 7 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m, n, p);
    }
    //! Returns the specified point from a 8 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    const double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m, n, p, q);
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
        return backingStore[static_cast<size_t>(arrayName)]->zIndexAndLayer(hIndex, layer);
    }

    //! Direct access top the underlying data array.
    const ModelArray& data() const { return *backingStore[static_cast<size_t>(arrayName)]; }
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
    const S& backingStore;
};

/*!
 * @brief A class which provides indirect access to ModelArray.
 *
 * @details Uses the ModelComponent registry of ModelArrays to access data from
 * other ModelComponents. The class provides indexing, casting and access (by
 * reference) to the pointed to ModelArray. Here the returned data are by
 * non-const references, used for SharedArray in read-write (RW) mode.
 *
 * @tparam sh The ModelComponent::SharedArray to be referenced.
 */
template <auto arrayName, typename S> class ModelArrayRef<arrayName, S, RW> {
public:
    ModelArrayRef(const S& backingStoreIn)
        : backingStore(backingStoreIn)
    {
    }

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
    double& operator[](const ModelArray::MultiDim& dims)
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator[](dims);
    }
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned.
     *
     * @param index The one dimensional index of the target point.
     */
    double& operator[](size_t index) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator[](index);
    }
    //! Returns the specified point from a 1 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i);
    }
    //! Returns the specified point from a 2 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j);
    }
    //! Returns the specified point from a 3 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k);
    }
    //! Returns the specified point from a 4 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k, size_t l) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l);
    }
    //! Returns the specified point from a 5 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m);
    }
    //! Returns the specified point from a 6 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m, n);
    }
    //! Returns the specified point from a 7 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m, n, p);
    }
    //! Returns the specified point from a 8 dimensional ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned.
    double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const
    {
        return backingStore[static_cast<size_t>(arrayName)]->operator()(i, j, k, l, m, n, p, q);
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
    double& zIndexAndLayer(size_t hIndex, size_t layer)
    {
        return backingStore[static_cast<size_t>(arrayName)]->zIndexAndLayer(hIndex, layer);
    }

    //! Direct access top the underlying data array.
    ModelArray& data() const { return *backingStore[static_cast<size_t>(arrayName)]; }
    //! Cast the reference class to a real reference to the referenced ModelArray.
    operator ModelArray&() const { return data(); }

private:
    const S& backingStore;
};
}
#endif /* MODELARRAYREF2_HPP */
