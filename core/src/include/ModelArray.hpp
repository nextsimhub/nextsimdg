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
#include <vector>

namespace Nextsim {

const static int DGdegree = 0; // TODO: Replace with the same source as the dynamics
const static int CellDoF = 1;
// TODO: (DGdegree == 0 ? 1 : (DGdegree == 1 ? 3 : (DGdegree == 2 ? 6 : -1)));
const static Eigen::StorageOptions majority = DGdegree == 0 ? Eigen::ColMajor : Eigen::RowMajor;

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
    enum class Type {
        H,
        U,
        V,
        Z,
        DG,
    };

    static const std::map<Type, std::string> typeNames;

    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, majority> DataType;
    //    typedef std::vector<double> DataType;

    class Component : Eigen::Array<double, 1, Eigen::Dynamic> {
    public:
        Component()
            : Eigen::Array<double, 1, Eigen::Dynamic>()
        {
        }
        Component(const Eigen::Array<double, 1, Eigen::Dynamic>& other)
            : Eigen::Array<double, 1, Eigen::Dynamic>(other)
        {
        }
    };

    static ModelArray HField(const std::string& name) { return ModelArray(Type::H, name); }
    static ModelArray UField(const std::string& name) { return ModelArray(Type::U, name); }
    static ModelArray VField(const std::string& name) { return ModelArray(Type::V, name); }
    static ModelArray ZField(const std::string& name) { return ModelArray(Type::Z, name); }

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
    ModelArray(const Type type, const std::string& name);
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

    typedef std::vector<size_t> Dimensions;

    //! Returns the number of dimensions of this type of ModelArray.
    size_t nDimensions() const { return nDimensions(type); }
    //! Returns the number of dimensions of the specified type of ModelArray.
    static size_t nDimensions(Type type) { return m_dims.at(type).size(); }
    //! Returns a vector<size_t> of the size of each dimension of this type of ModelArray.
    const Dimensions& dimensions() const { return dimensions(type); }
    //! Returns a vector<size_t> of the size of each dimension of the specified type of ModelArray.
    static const Dimensions& dimensions(Type type) { return m_dims.at(type); }
    //! Returns the total number of elements of this type of ModelArray.
    size_t size() const { return size(type); }
    //! Returns the total number of elements of the specified type of ModelArray.
    static size_t size(Type type) { return m_sz.at(type); }
    //! Returns the size of the data array of this object.
    size_t trueSize() const { return m_data.rows(); }

    //! Returns a read-only pointer to the underlying data buffer.
    const double* getData() const { return m_data.data(); }
    //! Retuns the (enum of) the ModelArray::Type of this.
    const Type getType() const { return type; }

    /*!
     * @brief Sets the number and size of the dimensions of a specified type of
     * ModelArray.
     *
     * @param type The type of array the dimensions are to be specified for.
     * @param dim The per-dimension size to be set.
     */
    static void setDimensions(Type, const Dimensions&);
    /*!
     * @brief Sets the number and size of the dimensions of this type of ModelArray.
     *
     * @details Sets the number and size of the dimensions of this type of
     * ModelArray. The data buffer of this object is then resized to match the
     * new definition.
     *
     * @param dim The per-dimension size to be set.
     */
    void setDimensions(const Dimensions& dims)
    {
        setDimensions(type, dims);
        resize();
    }

    //! Conditionally updates the size of the object data buffer to match the
    //! class specification.
    void resize()
    {
        if (size() != trueSize())
            m_data.resize(m_sz.at(type), Eigen::NoChange);
    }

    //! Returns the name of the object.
    const std::string& name() const { return m_name; }

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
    template <typename T, typename I>
    static inline T indexrHelper(const T* dims, T cdim, T& stride, I first)
    {
        stride *= dims[cdim];
        return first;
    }

    template <typename T, typename I, typename... Args>
    static inline T indexrHelper(const T* dims, T cdim, T& stride, I first, Args... args)
    {
        T lower = indexrHelper(dims, cdim + 1, stride, args...);
        T incr = stride * first;
        stride *= dims[cdim];
        return incr + lower;
    }

    template <typename T, typename I> static inline T indexr(const T* dims, I first)
    {
        return first;
    }

    template <typename T, typename I> static inline T indexr(const T* dims, I first, I second)
    {
        return first * dims[1] + second;
    }

    template <typename T, typename... Args> static inline T indexr(const T* dims, Args... args)
    {
        T stride = 1;
        return indexrHelper(dims, static_cast<T>(0), stride, args...);
    }

    static size_t indexr(const ModelArray::Dimensions& loc, const size_t* dims)
    {
        size_t loc8[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        for (char i = 0; i < loc.size(); ++i) {
            loc8[i] = loc[i];
        }

        return indexr(dims, loc8[0], loc8[1], loc8[2], loc8[3], loc8[4], loc8[5], loc8[6], loc8[7]);
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
    const double& operator[](const Dimensions& dims) const;

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
    double& operator[](const Dimensions&);
    //! Returns the specified point from a ModelArray. If the
    //! object holds discontinuous Galerkin components, only the cell averaged
    //! value is returned. Non-const version.
    template <typename... Args> double& operator()(Args... args)
    {
        return const_cast<double&>(std::as_const(*this)(args...));
    }

    /*!
     * @brief Accesses the full Discontinuous Galerkin coefficient vector at
     * the indexed location.
     *
     * @param i one-dimensional index of the target point.
     */
    Component components(size_t i) { return Component(m_data.col(i)); }

    /*!
     * @brief Accesses the full Discontinuous Galerkin coefficient vector at the specified location.
     *
     * @param dims indexing argument of the target point.
     */
    Component components(const Dimensions& loc);

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

    //! Returns the number of discontinuous Galerkin components held in this
    //! type of ModelArray.
    inline size_t nComponents() const { return nComponents(type); }
    //! Returns the number of discontinuous Galerkin components held in the
    //! specified type of ModelArray.
    inline static size_t nComponents(const Type type) { return (hasDoF(type)) ? CellDoF : 1; }
    //! Returns whether this type of ModelArray has additional discontinuous
    //! Galerkin components.
    inline bool hasDoF() const { return hasDoF(type); }
    //! Returns whether the specified type of ModelArray has additional
    //! discontinuous Galerkin components.
    inline static bool hasDoF(const Type type) { return type == Type::DG; }

private:
    class SizeMap {
    public:
        SizeMap()
            : m_sizes({ { Type::H, 0 }, { Type::U, 0 }, { Type::V, 0 }, { Type::Z, 0 } })
        {
        }
        size_t& at(const Type& type) { return m_sizes.at(type); }
        const size_t& at(const Type& type) const { return m_sizes.at(type); }

        size_t& operator[](const Type& type) { return m_sizes[type]; }
        size_t& operator[](Type&& type) { return m_sizes[type]; }

        size_t size() const noexcept { return m_sizes.size(); }

    private:
        std::map<Type, size_t> m_sizes;
    };
    static SizeMap m_sz;

    class DimensionMap {
    public:
        DimensionMap()
            : m_dimensions({ { Type::H, { 0 } }, { Type::U, { 0 } }, { Type::V, { 0 } },
                { Type::Z, { 0, 1 } } })
        {
        }
        Dimensions& at(const Type& type) { return m_dimensions.at(type); }
        const Dimensions& at(const Type& type) const { return m_dimensions.at(type); }

        Dimensions& operator[](const Type& type) { return m_dimensions[type]; }
        Dimensions& operator[](Type&& type) { return m_dimensions[type]; }

        size_t size() const noexcept { return m_dimensions.size(); }

    private:
        std::map<Type, Dimensions> m_dimensions;
    };
    static DimensionMap m_dims;
    DataType m_data;
    std::string m_name;
};

typedef ModelArray HField;
typedef ModelArray UField;
typedef ModelArray VField;
typedef ModelArray ZField;

// ModelArray arithmetic with doubles
ModelArray operator+(const double&, const ModelArray&);
ModelArray operator-(const double&, const ModelArray&);
ModelArray operator*(const double&, const ModelArray&);
ModelArray operator/(const double&, const ModelArray&);
} /* namespace Nextsim */

#endif /* MODELARRAY_HPP */
