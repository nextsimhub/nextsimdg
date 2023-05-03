/*!
 * @file ModelData.hpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAY_HPP
#define MODELARRAY_HPP

#include <Eigen/Core>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Nextsim {

/*
 * Set the storage order to row major. This matches with DGVector when there is
 * more than one DG component. If there is only one DG component (the finite
 * element component), then the order of the data in the buffer is the same,
 * and data can be safely transferred between the buffers underlying the two
 * types. The physics code will never directly touch the spatially varying
 * components, so the choice of storage order should not matter.
 */
const static Eigen::StorageOptions majority = Eigen::RowMajor;

/*!
 * @brief A class that holds the array data for the model.
 *
 * @details This class holds n-dimensional data for the model. It is designed
 * to be a light (not transparent) wrapper around Eigen::Array. The underlying
 * Array is two dimensional, with the first dimension acting like the one
 * dimensional vector backing an n-dimensional array type. All of the n
 * dimensional indexing is collapsed onto the first dimension. Any access to
 * components for Discontinuous Galerkin variables is done through the second
 * dimension. For finite volume variables, there is one component in the second
 * dimension.
 *
 * Sizes and dimensionality as calculated by the static and member functions
 * only pertain to the grid, with the DG components being ignored. For a
 * degree-2 DG variable the size in memory would be 6 times that reported by
 * the size() function.
 */
class ModelArray {
public:
    // Forward defines make Eclipse less red and squiggly
    enum class Type;
    enum class Dimension;

#include "include/ModelArrayDetails.hpp"

    static const int N_DEFINED_DIMENSIONS = static_cast<int>(Dimension::COUNT);

    struct DimensionSpec {
        std::string name;
        size_t length;
    };
    typedef std::map<Type, std::vector<Dimension>> TypeDimensions;

    //! The dimensions that make up each defined type. Defined in ModelArrayDetails.cpp
    static TypeDimensions typeDimensions;
    //! The name and length of each dimension that is defined
    static std::map<Dimension, DimensionSpec> definedDimensions;
    //! The name of each type of ModelArray
    static const std::map<Type, std::string> typeNames;
    // The dimension that defines the components of each ModelArray type, if any
    static const std::map<Type, Dimension> componentMap;

    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, majority> DataType;

    typedef DataType::RowXpr Component;
    typedef DataType::ConstRowXpr ConstComponent;

    /*!
     * Construct an unnamed ModelArray of Type::H
     */
    ModelArray();
    /*!
     * @brief Construct a ModelArray of the given type and name
     *
     * @param type The ModelArray::Type for the new object.
     * @param name The name of the new object.
     */
    ModelArray(const Type type);
    //! Copy constructor
    ModelArray(const ModelArray&);
    virtual ~ModelArray() {};

    //! Copy assignment operator
    ModelArray& operator=(const ModelArray&);
    /*!
     * @brief Assigns a double value to all elements of the object.
     *
     * @param val The value to be assigned.
     */
    ModelArray& operator=(const double& val);

    // ModelArray arithmetic
    //! In place addition of another ModelArray
    ModelArray& operator+=(const ModelArray& b)
    {
        m_data += b.m_data;
        return *this;
    }
    //! In place subtraction of another ModelArray
    ModelArray& operator-=(const ModelArray& b)
    {
        m_data -= b.m_data;
        return *this;
    }
    //! In place multiplication by another ModelArray
    ModelArray& operator*=(const ModelArray& b)
    {
        m_data *= b.m_data;
        return *this;
    }
    //! In place division by another ModelArray
    ModelArray& operator/=(const ModelArray& b)
    {
        m_data /= b.m_data;
        return *this;
    }

    //! In place addition of a double
    ModelArray& operator+=(double b)
    {
        m_data += b;
        return *this;
    }
    //! In place subtraction of a double
    ModelArray& operator-=(double b)
    {
        m_data -= b;
        return *this;
    }
    //! In place multiplication by a double
    ModelArray& operator*=(double b)
    {
        m_data *= b;
        return *this;
    }
    //! In place division by a double
    ModelArray& operator/=(double b)
    {
        m_data /= b;
        return *this;
    }

    //! Returns a ModelArray containing the per-element sum of the
    //! object and the provided ModelArray.
    ModelArray operator+(const ModelArray&) const;
    //! Returns a ModelArray containing the per-element difference between the
    //! object and the provided ModelArray.
    ModelArray operator-(const ModelArray&) const;
    //! Returns a ModelArray containing the per-element product of the
    //! object and the provided ModelArray.
    ModelArray operator*(const ModelArray&) const;
    //! Returns a ModelArray containing the per-element ratio between the
    //! object and the provided ModelArray.
    ModelArray operator/(const ModelArray&) const;
    // Returns a ModelArray containing the element-wise negation of this.
    ModelArray operator-() const;

    //! Returns a ModelArray with a constant added to every element of the object.
    ModelArray operator+(const double&) const;
    //! Returns a ModelArray with a constant subtracted from every element of the object.
    ModelArray operator-(const double&) const;
    //! Returns a ModelArray with every element of the object multiplied by a constant.
    ModelArray operator*(const double&) const;
    //! Returns a ModelArray with every element of the object divided by a constant.
    ModelArray operator/(const double&) const;

    /*!
     * @brief Calculates element-wise maximum of the data and the given scalar value.
     * @param max the maximum value of the resultant array.
     */
    ModelArray max(double max) const;
    /*!
     * @brief Calculates element-wise minimum of the data and the given scalar value.
     * @param min the minimum value of the resultant array.
     */
    ModelArray min(double min) const;
    /*!
     * @brief Calculates element-wise maximum of the data and the given second array.
     * @param maxArr the array containing the maximum values.
     */
    ModelArray max(const ModelArray& maxArr) const;
    /*!
     * @brief Calculates element-wise minimum of the data and the given second array.
     * @param minArr the array containing the minimum values.
     */
    ModelArray min(const ModelArray& minArr) const;

    /*!
     * @brief Clamps the values in the array to the given maximum.
     * @param max the maximum value of the final array.
     */
    ModelArray& clampAbove(double max);
    /*!
     * @brief Clamps the values in the array to the given minimum.
     * @param min the minimum value of the final array.
     */
    ModelArray& clampBelow(double min);
    /*!
     * @brief Clamps the values in the array to maximum values in the given array.
     * @param maxArr the array of clamp maximum target values.
     */
    ModelArray& clampAbove(const ModelArray& maxArr);
    /*!
     * @brief Clamps the values in the array to minimum values in the given array.
     * @param minArr the array of clamp minimum target values.
     */
    ModelArray& clampBelow(const ModelArray& minArr);

    typedef std::vector<size_t> MultiDim;

    //! Returns the number of dimensions of this type of ModelArray.
    size_t nDimensions() const { return nDimensions(type); }
    //! Returns the number of dimensions of the specified type of ModelArray.
    static size_t nDimensions(Type type) { return m_dims.at(type).size(); }
    //! Returns a vector<size_t> of the size of each dimension of this type of ModelArray.
    const MultiDim& dimensions() const { return dimensions(type); }
    //! Returns a vector<size_t> of the size of each dimension of the specified type of ModelArray.
    static const MultiDim& dimensions(Type type) { return m_dims.at(type); }
    //! Returns the total number of elements of this type of ModelArray.
    size_t size() const { return size(type); }
    //! Returns the total number of elements of the specified type of ModelArray.
    static size_t size(Type type) { return m_sz.at(type); }
    //! Returns the size of the data array of this object.
    size_t trueSize() const { return m_data.rows(); }

    //! Returns a read-only pointer to the underlying data buffer.
    const double* getData() const { return m_data.data(); }

    //! Returns a const reference to the Eigen data
    const DataType& data() const { return m_data; }
    //! Returns the (enum of) the ModelArray::Type of this.
    const Type getType() const { return type; }

    /*!
     * @brief Sets the number and size of the dimensions of a specified type of
     * ModelArray.
     *
     * @param type The type of array the dimensions are to be specified for.
     * @param dim The per-dimension size to be set.
     */
    static void setDimensions(Type, const MultiDim&);
    /*!
     * @brief Sets the number and size of the dimensions of this type of ModelArray.
     *
     * @details Sets the number and size of the dimensions of this type of
     * ModelArray. The data buffer of this object is then resized to match the
     * new definition.
     *
     * @param dim The per-dimension size to be set.
     */
    void setDimensions(const MultiDim& dims)
    {
        setDimensions(type, dims);
        resize();
    }

    /*!
     * @brief Sets the length of an individual dimension before propagating it
     * to the defined array types.
     *
     * @param dim The dimension to be altered.
     * @param length The new length of the dimension.
     */
    static void setDimension(Dimension dim, size_t length);

    //! Conditionally updates the size of the object data buffer to match the
    //! class specification.
    void resize()
    {
        if (size() != trueSize()) {
            if (hasDoF(type)) {
                m_data.resize(m_sz.at(type), definedDimensions.at(componentMap.at(type)).length);
            } else {
                m_data.resize(m_sz.at(type), Eigen::NoChange);
            }
        }
    }

    /*!
     * @brief Sets the value of every element in the object to the provided value.
     *
     * @param value The new value for every element.
     */
    void setData(double value);
    /*!
     * @brief Reads and sets data from a raw buffer of double data.
     *
     * @details The given pointer is the address of the first element read.
     * Subsequent doubles in memory are read until every element in the object
     * has been read. Ensure the provided buffer contains sufficient data to
     * fill every element of the object, as no further bounds checking is performed.
     *
     * @param pData The pointer to the data buffer to be read.
     */
    void setData(const double* pData);
    /*!
     * @brief Reads and sets data from an instance of the underlying data class
     * (Eigen::Array).
     *
     * @param data The data object to be copied from.
     */
    void setData(const DataType& data);
    /*!
     * @brief Reads and sets data from another ModelArray.
     *
     * @details This function reads the data buffer only, and performs no
     * bounds checking on the ModelArray argument. Unlike the assignment
     * operator or copy constructor, it does not change the size or type of the
     * object.
     *
     * @param source The object to be copied from.
     */
    void setData(const ModelArray& source);

private:
    // Fast special case for 1-d indexing
    template <typename T, typename I> static inline T indexr(const T* dims, I first)
    {
        return static_cast<T>(first);
    }

    // Fast special case for 2-d indexing
    template <typename T, typename I> static inline T indexr(const T* dims, I first, I second)
    {
        return first + second * dims[0];
    }

    // Indices as separate function parameters
    template <typename T, typename I, typename... Args>
    static inline T indexr(const T* dims, I first, Args... args)
    {
        std::initializer_list<I> loc { first, args... };
        return indexrHelper(dims, loc);
    }

    // Indices as a Dimensions object
    template <typename T> static T indexr(const T* dims, const ModelArray::MultiDim& loc)
    {
        return indexrHelper(dims, loc);
    }

    // Generic index generator that will work on any container
    template <typename T, typename C> static T indexrHelper(const T* dims, const C& loc)
    {
        size_t ndims = loc.size();
        T stride = 1;
        T ii = 0;
        auto iloc = begin(loc);
        for (size_t dim = 0; dim < ndims; ++dim) {
            ii += stride * (*iloc++);
            stride *= dims[dim];
        }
        return ii;
    }

public:
    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned. Here the data is a const reference.
     *
     * @param i The one dimensional index of the target point.
     */
    const double& operator[](size_t i) const { return m_data(i, 0); }
    /*!
     * @brief Returns the data at the indices.
     *
     * @details The argument is a list of dimension indices (actually a
     * std::vector<size_t>). The number of dimensions provided can be lower
     * than that of the ModelArray type. If the object holds discontinuous
     * Galerkin components, only the cell averaged value is returned. Here the
     * data is a const reference.
     *
     * @param dims The indices of the target point.
     */
    const double& operator[](const MultiDim& dims) const;

    /*!
     * @brief Returns the data at the given set of indices
     */
    template <typename... Args> const double& operator()(Args... args) const
    {
        return (*this)[indexr(dimensions().data(), args...)];
    }

    /*!
     * @brief Returns the data at the specified one dimensional index.
     *
     * @details The argument is used to directly index the data buffer. If the
     * object holds discontinuous Galerkin components, only the cell averaged
     * value is returned. Here the data is a reference and so can be used for
     * assignment.
     *
     * @param i The one dimensional index of the target point.
     */
    double& operator[](size_t i) { return const_cast<double&>(std::as_const(*this)(i)); }
    /*!
     * @brief Returns the data at the indices.
     *
     * @details The argument is a list of dimension indices (actually a
     * std::vector<size_t>). The number of dimensions provided can be lower
     * than that of the ModelArray type. If the object holds discontinuous
     * Galerkin components, only the cell averaged value is returned. Here the
     * data is a reference and so can be used for assignment.
     *
     * @param dims The indices of the target point.
     */
    double& operator[](const MultiDim&);
    //! Returns the specified point from a ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned. Non-const version.
    template <typename... Args> double& operator()(Args... args)
    {
        return const_cast<double&>(std::as_const(*this)(args...));
    }

    /*!
     * @brief Sets the number of components for DG & CG array types.
     *
     * @param cMap a map from ModelArray::Type to the number of components.
     */
    static void setNComponents(std::map<Type, size_t> cMap);

    /*!
     * @brief Sets the number of components for a single D/CG array type.
     *
     * @param type the DG or CG array type to set the number of components for.
     * @param nComp the number of components to be set.
     */
    static void setNComponents(Type type, size_t nComp)
    {
        if (hasDoF(type)) {
            definedDimensions.at(componentMap.at(type)).length = nComp;
        }
    }

    /*!
     * @brief Sets the number of components for this array's type.
     *
     * @param nComp the number of components to be set.
     */
    inline void setNComponents(size_t nComp) { setNComponents(type, nComp); }

    /*!
     * @brief Accesses the full Discontinuous Galerkin coefficient vector at
     * the indexed location.
     *
     * @param i one-dimensional index of the target point.
     */
    Component components(size_t i) { return m_data.row(i); }

    const ConstComponent components(size_t i) const { return m_data.row(i); }

    /*!
     * @brief Accesses the full Discontinuous Galerkin coefficient vector at the specified location.
     *
     * @param dims indexing argument of the target point.
     */
    Component components(const MultiDim& loc);
    const ConstComponent components(const MultiDim& loc) const;

    /*!
     * @brief Special access function for ZFields.
     *
     * @detail Index a ZField using an index from an HField of the same
     * horizontal extent and a layer index for the final dimension.
     *
     * @param hIndex the equivalent positional index in an HField array
     * @param layer the vertical layer to be accessed
     */
    double& zIndexAndLayer(size_t hIndex, size_t layer)
    {
        return this->operator[](zLayerIndex(hIndex, layer));
    }
    /*!
     * @brief Special access function for ZFields, const version.
     *
     * @detail Index a ZField using an index from an HField of the same
     * horizontal extent and a layer index for the final dimension.
     *
     * @param hIndex the equivalent positional index in an HField array
     * @param layer the vertical layer to be accessed
     */
    const double& zIndexAndLayer(size_t hIndex, size_t layer) const
    {
        return this->operator[](zLayerIndex(hIndex, layer));
    }

    /*!
     * @brief Returns the index for a given set of multi-dimensional location for this array's type.
     *
     * @param loc The multi-dimensional location to return the index for.
     */
    size_t indexFromLocation(const MultiDim& loc) const { return indexFromLocation(type, loc); }

    /*!
     * @brief Returns the index for a given set of multi-dimensional location for the given Type.
     *
     * @param type The type to act on.
     * @param loc The multi-dimensional location to return the index for.
     */
    static size_t indexFromLocation(Type type, const MultiDim& loc);

    /*!
     * @brief Returns the multi-dimensional location for a given index for this array's type.
     *
     * @param index The index to return the multi-dimensional location for.
     */
    MultiDim locationFromIndex(size_t index) const { return locationFromIndex(type, index); }

    /*!
     * @brief Returns the multi-dimensional location for a given index for the given Type.
     *
     * @param type The type to act on.
     * @param index The index to return the multi-dimensional location for.
     */
    static MultiDim locationFromIndex(Type type, size_t index);

    //! Returns the number of discontinuous Galerkin components held in this
    //! type of ModelArray.
    inline size_t nComponents() const { return nComponents(type); }
    //! Returns the number of discontinuous Galerkin components held in the
    //! specified type of ModelArray.
    inline static size_t nComponents(const Type type)
    {
        return (hasDoF(type)) ? definedDimensions.at(componentMap.at(type)).length : 1;
    }
    //! Returns whether this type of ModelArray has additional discontinuous
    //! Galerkin components.
    inline bool hasDoF() const { return hasDoF(type); }
    //! Returns whether the specified type of ModelArray has additional
    //! discontinuous Galerkin components.
    static bool hasDoF(const Type type);

protected:
    Type type;

    /*!
     * @brief Special access function for ZFields, common implementation version.
     *
     * @detail Index a ZField using an index from an HField of the same
     * horizontal extent and a layer index for the final dimension.
     *
     * @param hIndex the equivalent positional index in an HField array
     * @param layer the vertical layer to be accessed
     */
    size_t zLayerIndex(size_t hIndex, size_t layer) const
    {
        return hIndex * dimensions()[nDimensions() - 1] + layer;
    }

private:
    static bool areMapsInvalid;
    static void validateMaps();
    class SizeMap {
    public:
        SizeMap();
        size_t& at(const Type& type) { return m_sizes.at(type); }
        const size_t& at(const Type& type) const { return m_sizes.at(type); }

        size_t& operator[](const Type& type) { return m_sizes[type]; }
        size_t& operator[](Type&& type) { return m_sizes[type]; }

        size_t size() const noexcept { return m_sizes.size(); }

        void validate();

    protected:
        std::map<Type, size_t> m_sizes;
    };
    static SizeMap m_sz;

    class DimensionMap {
    public:
        DimensionMap();
        MultiDim& at(const Type& type) { return m_dimensions.at(type); }
        const MultiDim& at(const Type& type) const { return m_dimensions.at(type); }

        MultiDim& operator[](const Type& type) { return m_dimensions[type]; }
        MultiDim& operator[](Type&& type) { return m_dimensions[type]; }

        size_t size() const noexcept { return m_dimensions.size(); }

        void validate();

    private:
        std::map<Type, MultiDim> m_dimensions;
    };
    static DimensionMap m_dims;
    DataType m_data;
};

#include "include/ModelArrayTypedefs.hpp"

// ModelArray arithmetic with doubles
ModelArray operator+(const double&, const ModelArray&);
ModelArray operator-(const double&, const ModelArray&);
ModelArray operator*(const double&, const ModelArray&);
ModelArray operator/(const double&, const ModelArray&);
} /* namespace Nextsim */

#endif /* MODELARRAY_HPP */
