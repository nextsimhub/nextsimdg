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

class DevGridIO;

class DevGrid : public IStructure {
public:
    DevGrid()
    : cursor(*this)
    , pio(nullptr)
    {
        processedStructureName = ourStructureName;
    }

    virtual ~DevGrid()
    {
        if (pio) {
            delete pio;
        }
    }

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

    class IDevGridIO {
    public:
        IDevGridIO(DevGrid& grid)
            : grid(&grid)
        {}
        virtual ~IDevGridIO() = default;
        virtual void init(std::vector<ElementData>& dg, const std::string& filePath) const = 0;
        virtual void dump(const std::vector<ElementData>& dg, const std::string& fielPath) const = 0;
    protected:
        DevGrid* grid;
    };

    IDevGridIO* pio;

    void setIO(IDevGridIO* p) { pio = p; }
private:
    const static std::string ourStructureName;
    const static std::string xDimName;
    const static std::string yDimName;

    std::vector<ElementData> data;

    std::vector<ElementData>::iterator iCursor;

    friend DevGridIO;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVGRID_HPP */
