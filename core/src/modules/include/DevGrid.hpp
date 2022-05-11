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
#include "include/IDevGridIO.hpp"
#include "include/ModelState.hpp"
#include "include/PrognosticElementData.hpp"

#include <map>

namespace Nextsim {

class DevGridIO;

//! A class to hold a grid of ElementData instances in a fixed sized square grid.
class DevGrid : public IStructure {
public:
    DevGrid()
        : pio(nullptr)
    {
    }

    //! Destructor. The lifetime of pio should be the lifetime of the instance.
    virtual ~DevGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    const static int nx;
    const static std::string structureName;

    // Read/write override functions
    void init(const std::string& filePath) override;

    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }

    void dump(const std::string& filePath) const override;

    void dumpModelState(const ModelState& state, const std::string& filePath) const override
    {
        if (pio)
            pio->dumpModelState(state, filePath);
    }
    std::string structureType() const override { return structureName; };

    int nIceLayers() const override { return 1; };

    // Cursor manipulation override functions
    int resetCursor() override;
    bool validCursor() const override;
    ElementData& cursorData() override;
    const ElementData& cursorData() const override;
    void incrCursor() override;

    //! Sets the pointer to the class that will perform the IO. Should be an instance of DevGridIO
    void setIO(IDevGridIO* p) { pio = p; }

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

private:
    std::vector<ElementData> data;

    std::vector<ElementData>::iterator iCursor;

    IDevGridIO* pio;

    friend DevGridIO;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVGRID_HPP */
