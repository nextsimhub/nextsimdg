/*!
 * @file ModelData.hpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELARRAY_HPP
#define CORE_SRC_INCLUDE_MODELARRAY_HPP

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace Nextsim {

class ModelArray {
public:
    enum class Type {
        H,
        U,
        V,
        Z,
    };

    static const std::map<Type, std::string> typeNames;

    static ModelArray HField(const std::string& name) { return ModelArray(Type::H, name); }
    static ModelArray UField(const std::string& name) { return ModelArray(Type::U, name); }
    static ModelArray VField(const std::string& name) { return ModelArray(Type::V, name); }
    static ModelArray ZField(const std::string& name) { return ModelArray(Type::Z, name); }

    ModelArray();
    ModelArray(const ModelArray&);
    virtual ~ModelArray() {};

    ModelArray& operator=(const ModelArray&);

    typedef std::vector<size_t> Dimensions;

    size_t nDimensions() const { return nDimensions(type); }
    static size_t nDimensions(Type type) { return m_dims.at(type).size(); }
    const Dimensions& dimensions() const { return dimensions(type); }
    static const Dimensions& dimensions(Type type) { return m_dims.at(type); }
    size_t size() const { return size(type); }
    size_t trueSize() const { return m_data.size(); }
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
            m_data.resize(m_sz.at(type));
    }

    const std::string& name() const { return m_name; }

    void setData(const double* pData);
    void setData(const std::vector<double>&);

    std::vector<double>::iterator begin() { return m_data.begin(); }
    std::vector<double>::iterator end() { return m_data.end(); }

    const double& operator[](size_t i) const { return m_data[i]; }
    const double& operator[](const Dimensions& dims) const;

    const double& operator()(size_t i) const { return m_data[i]; }
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
    double& operator[](size_t i) { return m_data[i]; }
    double& operator[](const Dimensions&);
    //! Multi (one) dimensional indexing
    double& operator()(size_t i) { return m_data[i]; }
    double& operator()(size_t i, size_t j);
    double& operator()(size_t i, size_t j, size_t k);
    double& operator()(size_t i, size_t j, size_t k, size_t l);
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);
    double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p);
    double& operator()(
        size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p, size_t q);

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
    std::vector<double> m_data;
    std::string m_name;
};

typedef ModelArray HField;
typedef ModelArray UField;
typedef ModelArray VField;
typedef ModelArray ZField;

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELARRAY_HPP */
