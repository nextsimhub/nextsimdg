/*!
 * @file DevGrid.hpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DEVGRID_HPP
#define CORE_SRC_INCLUDE_DEVGRID_HPP

#include "include/IStructure.hpp"

#include "include/ElementData.hpp"
#include "include/PrognosticData.hpp"

#include <map>

namespace Nextsim {

class DevGrid : public IStructure {
public:
    DevGrid()
    : cursor(*this)
    {
        processedStructureName = ourStructureName;
    }

    virtual ~DevGrid() = default;

    const static int nx;

    // Read/write override functions
    void init(const std::string& filePath) override;
    void dump(const std::string& filePath) const override;

    // Cursor manipulation override functions
    int resetCursor() override;
    bool validCursor() const override;
    ElementData& cursorData() override;
    const ElementData& cursorData() const override;
    void incrCursor() override;

    class Cursor : public IStructure::Cursor {
    public:
        Cursor(DevGrid& dg)
            : owner(dg)
        {
        }

        ~Cursor() = default;
        IStructure& operator=(const int) const override;
        operator bool() const override;
        ElementData& operator*() const override;
        ElementData* operator->() const override;
        IStructure& operator++() const override;

    private:
        DevGrid& owner;
    };

    const Cursor cursor;

protected:
    /*!
     * @brief Initializes the structure based on the contents of the structure
     * group of the input file.
     *
     * @param grp The NetCDF group instance that holds the structure
     * information.
     */
    void init(const netCDF::NcGroup& grp);
    /*!
     * @brief Initializes the structure of the IStructure from metadata.
     *
     * @param metaGroup The NetCDF group instance holding the structure
     * metadata.
     */
    void initMeta(const netCDF::NcGroup& metaGroup);
    /*!
     * @brief Initializes the contents of the IStructure from data.
     *
     * @param dataGroup The NetCDF group instance holding the data.
     */
    void initData(const netCDF::NcGroup& dataGroup);
    /*!
     * @brief Dumps the data and metadata to two sub-groups.
     *
     * @details The structure metadata will be dumped to immediate subnodes
     * with names specified by the ::metadataNodeName and ::dataNodeName public
     * member variables.
     *
     * @param headGroup The top-level node to hold the metadata and data nodes.
     */
    void dump(netCDF::NcGroup& headGroup) const;
    /*!
     * @brief Dumps the structural metadata to a netCDF node.
     *
     * @param metaGroup The top-level node to write the metadata to.
     */
    void dumpMeta(netCDF::NcGroup& metaGroup) const;
    /*!
     * @brief Dumps the prognostic data to a netCDF node.
     *
     * @param dataGroup The top-level node to write the data to.
     */
    void dumpData(netCDF::NcGroup& dataGroup) const;

private:
    const static std::string ourStructureName;
    const static std::string xDimName;
    const static std::string yDimName;

    // pointer to a member of PrognosticData that takes no arguments and returns a
    // double. See https://isocpp.org/wiki/faq/pointers-to-members#typedef-for-ptr-to-memfn
    typedef double (PrognosticData::*ProgDoubleFn)() const;

    // Map between variable names and retrieval functions
    static const std::map<std::string, ProgDoubleFn> variableFunctions;

    std::vector<ElementData> data;

    std::vector<ElementData>::iterator iCursor;

    std::vector<double> gather(ProgDoubleFn pFunc) const;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVGRID_HPP */
