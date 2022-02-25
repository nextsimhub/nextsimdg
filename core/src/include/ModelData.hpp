/*!
 * @file ModelData.hpp
 *
 * @date Feb 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_MODELDATA_HPP
#define CORE_SRC_INCLUDE_MODELDATA_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace Nextsim {

class ModelData {
public:
    ModelData();
    ModelData(const std::string&);
    ModelData(const ModelData&);
    virtual ~ModelData() {};

    ModelData& operator=(const ModelData&);

    typedef std::vector<size_t> Dimensions;

    static size_t nDimensions() { return m_dims.size(); }
    static const Dimensions& dimensions() { return m_dims; }
    static size_t size() { return m_sz; }

    static void setDimensions(const Dimensions& dims);

    const std::string& name() const { return m_name; }

    void setData(const double* pData);
    void setData(const std::vector<double>&);

    std::vector<double>::iterator begin() { return m_data.begin(); }
    std::vector<double>::iterator end() { return m_data.end(); }

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

private:
    static size_t m_sz;
    static Dimensions m_dims;
    std::vector<double> m_data;
    std::string m_name;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_MODELDATA_HPP */
