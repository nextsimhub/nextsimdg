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
    DevGrid();
    virtual ~DevGrid();

    const static int nx;

    // Read/write override functions
    void initMeta(const netCDF::NcGroup& metaGroup) override;
    void initData(const netCDF::NcGroup& dataGroup) override;
    void dumpMeta(netCDF::NcGroup& metaGroup) const override;
    void dumpData(netCDF::NcGroup& dataGroup) const override;

    // Cursor manipulation override functions
    int resetCursor() override;
    bool validCursor() const override;
    ElementData& cursorData() override;
    const ElementData& cursorData() const override;
    void incrCursor() override;

    class Cursor : public IStructure::Cursor {
    public:
        Cursor(DevGrid&);
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
