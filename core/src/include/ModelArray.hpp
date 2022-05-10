/*!
 * @file ModelData.hpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELARRAY_HPP
#define CORE_SRC_INCLUDE_MODELARRAY_HPP

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

    ModelArray();
    ModelArray(const ModelArray&);
    virtual ~ModelArray() {};

    ModelArray& operator=(const ModelArray&);

    // ModelArray arithmetic
    ModelArray operator+(const ModelArray&) const;
    ModelArray operator-(const ModelArray&) const;
    ModelArray operator*(const ModelArray&) const;
    ModelArray operator/(const ModelArray&) const;

    typedef std::vector<size_t> Dimensions;

    //! Returns the number of dimensions of the physical grid
    size_t nDimensions() const { return nDimensions(type); }
    static size_t nDimensions(Type type) { return m_dims.at(type).size(); }
    const Dimensions& dimensions() const { return dimensions(type); }
    static const Dimensions& dimensions(Type type) { return m_dims.at(type); }
    size_t size() const { return size(type); }
    size_t trueSize() const { return m_data.rows(); }
    static size_t size(Type type) { return m_sz.at(type); }

    static void setDimensions(Type, const Dimensions&);
    void setDimensions(const Dimensions& dims)
    {
        setDimensions(type, dims);
        resize();
    }

    void resize()
    {
        if (size() != trueSize())
            m_data.resize(m_sz.at(type), Eigen::NoChange);
    }

    const std::string& name() const { return m_name; }

    void setData(double value);
    void setData(const double* pData);
    void setData(const DataType&);

    const double& operator[](size_t i) const { return m_data(i, 0); }
    const double& operator[](const Dimensions& dims) const;

    const double& operator()(size_t i) const { return (*this)[i]; }
    const double& operator()(size_t i, size_t j) const;
    const double& operator()(size_t i, size_t j, size_t k) const;
    const double& operator()(size_t i, size_t j, size_t k, size_t l) const;
    const double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;
    const double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;
    const double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p) const;
    const double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q) const;

    //! One dimensional indexing
    double& operator[](size_t i) { return const_cast<double&>(std::as_const(*this)(i)); }
    double& operator[](const Dimensions&);
    //! Multi (one) dimensional indexing
    double& operator()(size_t i) { return (*this)[i]; }
    double& operator()(size_t i, size_t j);
    double& operator()(size_t i, size_t j, size_t k);
    double& operator()(size_t i, size_t j, size_t k, size_t l);
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p);
    double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q);

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
    ModelArray(const Type, const std::string&);

    size_t zLayerIndex(size_t hIndex, size_t layer) const
    {
        return hIndex * dimensions()[nDimensions() - 1] + layer;
    }

    inline size_t nComponents() const { return nComponents(type); }
    inline static size_t nComponents(const Type type) { return (hasDoF(type)) ? CellDoF : 1; }
    inline bool hasDoF() const { return hasDoF(type); }
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

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELARRAY_HPP */
